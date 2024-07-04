#![allow(dead_code)]
#![allow(unused_variables)]
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]
#![allow(unused_assignments)]

use std::f32::consts::PI;

use crate::twoline::{Regime, TwoPhaseLine};

const G: f64 = 9.81; // gravity accelerator [,/s^2]
const GC: f64 = 9.8; // gravity constant [kg-m/kgf-s^2]

pub struct Horizontal {
    // process data
    pub WL: f64,     // liquid mass flow rate [kg/hr]
    pub WG: f64,     // Vapor mass flow rate [kg/hr]
    pub LoL: f64,    // Liquid density [kg/m^3]
    pub LoG: f64,    // Vapor density [kg/m^3]
    pub muL: f64,    // Liquid viscosity [cP] -> [kg/m-s]
    pub muG: f64,    // Vapor viscosity [cP] -> [kg/m-s]
    pub ST: f64,     // Liquid surface tension [dyne/cm] -> [kg/s^2]
    pub rough: f64,  // pipe absolute roughness [mm] -> [m]
    pub SF: f64,     // Safety factor [-]
    pub ID: f64,     // pipe inside diameter [in] -> [m]
    pub degree: f64, // degree,  Horizontal = 0, -Up / +Down

    // result
    pub regime_enum: Regime, // identify the flow regime(enum)
    pub flow_regime: String, // identify the flow regime(String)

    // state variable
    is_unit_change: bool, // is call the unit_transfer and transfer to unit
}

impl Horizontal {
    pub fn new(
        WL: f64,
        WG: f64,
        LoL: f64,
        LoG: f64,
        muL: f64,
        muG: f64,
        ST: f64,
        rough: f64,
        SF: f64,
        ID: f64,
        degree: f64,
    ) -> Self {
        Horizontal {
            WL,
            WG,
            LoL,
            LoG,
            muL,
            muG,
            ST,
            rough,
            SF,
            ID,
            degree,
            is_unit_change: false,
            regime_enum: Regime::NONE,
            flow_regime: String::from(""),
        }
    }
}

impl Horizontal {
    fn fhLL(&self, DD: f64, h: f64, XS: f64) -> f64 {
        let term1 = (2.0 * h - 1.0).acos(); // Eq. (13)
        let term2 = (1.0 - (2.0 * h - 1.0).powi(2)).sqrt(); // Eq. (14)
        let ALB = 0.25 * (std::f64::consts::PI - term1 + (2.0 * h - 1.0) * term2); // Eq. (10)
        let AGB = 0.25 * (term1 - (2.0 * h - 1.0) * term2); // Eq. (11)
        let SLB = std::f64::consts::PI - term1; // Eq. (12)
        let SGB = term1; // Eq. (13)
        let SiB = term2; // Eq. (14)
        let AB = std::f64::consts::PI / 4.0; // 相對於面積參考量 D^2 的無因次管截面積
        let ULB = AB / ALB; // Eq. (15)
        let UGB = AB / AGB; // Eq. (16)
        let DLB = 4.0 * ALB / SLB; // Eq. (6)
        let DGB = 4.0 * AGB / (SGB + SiB); // Eq. (6)
        let term3 = (self.LoL - self.LoG) * G * self.degree.sin(); // Eq. (9) Numerator
        let CG = 0.046; // 氣體摩擦因子關聯式中的常數 for turbulent flow
        let m = 0.2; // Eq. (5) 中的次幕 for turbulent flow
        let nuG = self.muG / self.LoG; // Dynamic Viscosity of Gas [Stoke]
        let area = std::f64::consts::PI * DD * DD / 4.0; // pipe area [m^2]
        let UX = self.WG / self.LoG / area / 3600.0; // Vapor Velocity [m/s]
        let term4 = 4.0 * CG / DD * (UX * DD / nuG).powf(-m) * (self.LoG * UX.powi(2) / 2.0); // Eq. (9) denominator
        let Y = term3 / term4; // 流體在流動方向之重力與壓降比值 (向上流動取負，向下流動取正，水平管 Y =0;

        let n = 0.2; // Eq. (5) 中的次幕 for turbulent flow
        let term5 = (ULB * DLB).powf(-n) * ULB.powi(2) * SLB / ALB;
        let term6 = (UGB * DGB).powf(-m) * UGB.powi(2) * (SGB / AGB + SiB / ALB + SiB / AGB);

        if XS == 0.0 {
            (4.0 * Y + term6 / term5).sqrt() // Eq. (7)
        } else {
            XS.powi(2) * term5 - term6 - 4.0 * Y
        }
    }

    fn ratioB(&self, X: f64, Y: f64) -> f64 {
        let hL = 0.5; // Kellogg modify to 0.35, but original still take 0.5
        let AB = std::f64::consts::PI / 4.0;
        let SGB = (2.0f64 * hL - 1.0f64).acos();
        let SLB = std::f64::consts::PI - (2.0 * hL - 1.0).acos();
        let SiB = (1.0 - (2.0 * hL - 1.0).powf(2.0)).sqrt();
        let ALB = 0.25 * (std::f64::consts::PI - SGB + (2.0 * hL - 1.0) * SiB);
        let AGB = 0.25 * (SGB - (2.0 * hL - 1.0) * SiB);
        let ULB = AB / ALB;
        let UGB = AB / AGB;
        let n = 0.2;
        let m = 0.2;
        let DLB = 4.0 * ALB / SLB;
        let DGB = 4.0 * AGB / (SGB + SiB);
        let term1 = (ULB * DLB).powf(-n) * ULB.powf(2.0) * SLB / ALB;
        let term2 = (UGB * DGB).powf(-m) * UGB.powf(2.0) * (SGB / AGB + SiB / ALB + SiB / AGB);
        let result = (X * X * term1 / (4.0 * Y + term2)).sqrt();
        result
    }
}

impl TwoPhaseLine for Horizontal {
    fn unit_transfer(&mut self) {
        if !self.is_unit_change {
            self.muL = self.muL * 0.001; // [cP] -> [kg/m-s]
            self.muG = self.muG * 0.001; // [cP] -> [kg/m-s]
            self.ID = self.ID * 2.54 / 100.0; // [in] -> [m]
            self.ST = self.ST * 1.019716213E-4; // [dyne/cm] -> [kgf/s^2]
            self.degree = self.degree * (PI as f64) / 180.0; // [degree] -> [rad]
            self.rough = self.rough / 1000.0; // [mm] -> [m]

            self.is_unit_change = true;
        }
    }

    fn flow_regime(&mut self) {
        // assume turbulent flow Eq.(8), see ref. 01
        let X = (self.WL / self.WG).powf(0.9)
            * (self.LoG / self.LoL).sqrt()
            * (self.muL / self.muG).powf(0.1);
        let mut hLa = 0.001f64; // hla: 波浪的平衡液位高 left initial value
        let mut hLb = 0.999f64; // hlb: 波浪的平衡液位高 right initial value
        let mut hLm; // hlm: hL mean value
        let eps = 1e-4; // eps: allowable tolerance

        // solve non-linear equation by Bisection Method
        loop {
            hLm = (hLa + hLb) / 2.0;
            if self.fhLL(self.ID, hLa, X) * self.fhLL(self.ID, hLm, X) < 0.0 {
                hLb = hLm;
            } else {
                hLa = hLm;
            }
            if (hLb - hLa).abs() <= eps {
                break;
            }
        }

        let hL = hLm; // hL 波浪的平衡液位高 [m]
        let area = std::f64::consts::PI * self.ID * self.ID / 4.0; // pipe area [m^2]
        let UX = self.WG / self.LoG / area / 3600.0; // Vapor Velocity [m/s]
        let UY = self.WL / self.LoL / area / 3600.0; // Liquid Velocity [m/s]
        let F = (self.LoG / (self.LoL - self.LoG)).sqrt() * UX
            / (self.ID * G * self.degree.cos()).sqrt(); // Froude Number Eq.(26)
        let C2 = 1.0 - hL; // Eq. (24)
        let term1 = (2.0 * hL - 1.0).acos(); // Eq. (13)
        let term2 = (1.0 - (2.0 * hL - 1.0).powi(2)).sqrt(); // Eq. (14)
        let dALB = term2;
        let ALB = 0.25 * (std::f64::consts::PI - term1 + (2.0 * hL - 1.0) * term2); // Eq. (10)
        let AGB = 0.25 * (term1 - (2.0 * hL - 1.0) * term2); // Eq. (11)
        let AB = std::f64::consts::PI / 4.0; // 相對於面積參考量 D^2 的無因次管截面積
        let ULB = AB / ALB; // Eq. (15)
        let UGB = AB / AGB; // Eq. (16)
        let ratio_a = ((F.powi(2) / C2.powi(2)) * UGB.powi(2) * dALB / AGB).sqrt();
        // Eq. (25) for Curve A

        // ratio C here
        let nuL = self.muL / self.LoL; // Dynamic Viscosity of Liquid [Stoke]
        let ReLS = self.ID * UY / nuL; // Liquid Slug Reynold Number [-]
        let K = F * ReLS.sqrt(); // Wavy flow dimensionless parameter [-]
        let S: f64 = 0.01; // 隱藏參數 [-]
        let ratio_c = K * ULB.sqrt() * UGB * S.sqrt() / 2.0; // Eq. (30) for Curve C

        // ratio B here
        let term3 = (self.LoL - self.LoG) * G * (self.degree).sin(); // Eq. (9) Numerator
        let CG = 0.046; // 氣體摩擦因子關聯式中的常數 for turbulent flow
        let CL = 0.046; // 氣體摩擦因子關聯式中的常數 for turbulent flow
        let n = 0.2; // Eq. (5) 中的次幕 for turbulent flow
        let m = 0.2; // Eq. (5) 中的次幕 for turbulent flow
        let nuG = self.muG / self.LoG; // Dynamic Viscosity of Gas [Stoke]
        let term4 =
            4.0 * CG / self.ID * (UX * self.ID / nuG).powf(-m) * (self.LoG * UX.powf(2.0) / 2.0); // Eq. (9) denominator
        let Y = term3 / term4; // 流體在流動方向之重力與壓降比值 (向上流動取負，向下流動取正，水平管 Y =0;
        let ratio_b = self.ratioB(X, Y);

        // ratio D here
        let t1 = 4.0 * CL / self.ID * (UY * self.ID / nuG).powf(-n) * self.LoL * UY.powf(2.0) / 2.0;
        let t2 = (self.LoL - self.LoG) * G * (self.degree).cos();
        let T2 = t1 / t2; // Eq. (37)

        let SiB = (1.0 - (2.0 * hL - 1.0).powf(2.0)).sqrt(); // Eq. (14)
        let SLB = std::f64::consts::PI - (2.0 * hL - 1.0).acos(); // Eq. (12)
        let DLB = 4.0 * ALB / SLB; // Eq. (6)
        let ratio_d = t2 * SiB * ULB.powf(2.0) * (ULB * DLB).powf(-n) / 8.0 / AGB;
        // Eq. (36)

        // EE here
        let UGScal = ((UY + G * (self.LoL - self.LoG) * self.ST / self.LoL.powf(2.0)).powf(0.25))
            * 1.15
            / 3.0;
        let EE = UX / UGScal;

        // judge regime by ratio
        if ratio_a <= 1.0 {
            // left side
            if ratio_c <= 1.0 {
                // down side
                self.regime_enum =
                    Regime::HorizontalStratifiedSmoothFlow(String::from("Stratified Smooth Flow"));
            } else {
                // top side
                self.regime_enum =
                    Regime::HorizontalStratifiedWavyFlow(String::from("Stratified Wavy Flow"));
            }
        } else {
            // right side
            if ratio_b <= 1.0 {
                self.regime_enum =
                    Regime::HorizontalAnnularDispersedFlow(String::from("Annular-Dispersed Flow"));
            } else {
                if ratio_d <= 1.0 {
                    if EE <= 1.0 {
                        self.regime_enum = Regime::HorizontalElongatedBubbleFlow(String::from(
                            "Elongated Bubble Flow",
                        ));
                    } else {
                        self.regime_enum = Regime::HorizontalIntermittentSlugFlow(String::from(
                            "Intermittent-Slug Flow",
                        ));
                    }
                } else {
                    self.regime_enum = Regime::HorizontalDispersedBubbleFlow(String::from(
                        "Dispersed Bubble Flow",
                    ));
                }
            }
        }
        self.flow_regime = match &self.regime_enum {
            Regime::HorizontalStratifiedSmoothFlow(v) => v.to_string(),
            Regime::HorizontalStratifiedWavyFlow(v) => v.to_string(),
            Regime::HorizontalAnnularDispersedFlow(v) => v.to_string(),
            Regime::HorizontalElongatedBubbleFlow(v) => v.to_string(),
            Regime::HorizontalIntermittentSlugFlow(v) => v.to_string(),
            Regime::HorizontalDispersedBubbleFlow(v) => v.to_string(),
            _ => String::from(""),
        };
    }

    fn model_cal(&mut self) {
        todo!()
    }
}
