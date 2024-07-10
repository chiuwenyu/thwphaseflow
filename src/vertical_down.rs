#![allow(dead_code)]
#![allow(unused_variables)]
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]
#![allow(unused_assignments)]

use std::f32::consts::PI;

use twoline::TwoPhaseLine;

use crate::twoline;
use crate::twoline::Regime;

const G: f64 = 9.81; // gravity accelerator [,/s^2]
const GC: f64 = 9.8; // gravity constant [kg-m/kgf-s^2]

pub struct VerticalDown {
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

    // for Annular Model
    pub LoTP: f64,  // Two-phase density [kg/m^3]
    pub UTP: f64,   // Two-phase velocity [m/s]
    pub alfaL: f64, // Liquid volume fraction
    pub Head: f64,  // 1.0 velocity head [kgf/cm^2]
    pub Pfric: f64, // Frictional pressure loss (kgf/cm^2/100m]
    pub Pgrav: f64, // Elevation Head Loss [kgf/cm^2/100m]
    pub Ef: f64,    // Erosion Factor [-]

    // for Bubble Model
    // LoTP, UTP, Head, Pfric, Pgrav, Ef same as Annular Model
    pub HL: f64, // Liquid Volume Fraction [-]

    // for Slug Model
    pub Loip: f64, // Two-phase density [kg/m^3]
    pub LoLS: f64, // Liquid Slug Density [kg/m^3]
    // HL, Head, Pfric, Pgrav, EF same as Annular Model

    // state variable
    is_unit_change: bool, // is call the unit_transfer and transfer to unit
}

impl crate::vertical_down::VerticalDown {
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
        crate::vertical_down::VerticalDown {
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
            LoTP: 0.0,
            UTP: 0.0,
            alfaL: 0.0,
            Head: 0.0,
            Pfric: 0.0,
            Pgrav: 0.0,
            Ef: 0.0,
            HL: 0.0,
            Loip: 0.0,
            LoLS: 0.0,
        }
    }

    fn get_uyc_from_curve_c(&self, x: f64) -> f64 {
        let area = std::f64::consts::PI * self.ID * self.ID / 4.0; // pipe area [m^2]
        let UL = self.WL / self.LoL / area / 3600.0; // Liquid Velocity [m/s]

        let term1 = 2.0 * (0.4 * self.ST / (self.LoL - self.LoG) / G).sqrt();
        let term2 = (self.LoL / self.ST).powf(0.6);
        let CL = 0.046;
        let n = 0.2; // Friction Factor parameter by Eq (11)
        let nuL = self.muL / self.LoL; // Liquid Kinetic Viscosity [m^2/s]
        let term3 = (2.0 / self.ID * CL * (self.ID / nuL).powf(-n)).powf(0.4);
        let term_b = term1 * term2 * term3;

        // Iterative method to find UM
        let mut UMi = UL + x; // initial value
        let trials = 100; // trial number
        let eps = 1e-6; // allowable tolerance
        let mut UMcal = 0.0; // reset calc value

        let mut i = 0;
        while i < trials {
            i += 1;
            let rt = 0.725 + 4.15 * (x / UMi).sqrt();
            let power = 2.0 * (3.0 - n) / 5.0;
            UMcal = (rt / term_b).powf(1.0 / power);
            let delta = (UMcal - UMi).abs();
            if delta > eps {
                UMi = UMcal;
            } else {
                break;
            }
        }

        let ULScal; // Curve B calculated ULS value
        if i < trials {
            ULScal = UMcal - x; // Convergence
        } else {
            ULScal = 0.0; // Divergence
        }

        ULScal
    }

    fn fanning(&self, Re: f64) -> f64 {
        // by Chen (1979)
        if Re < 2100.0 {
            16.0 / Re
        } else {
            let a = (self.rough / self.ID).powf(1.1098) / 2.8257 + (7.149 / Re).powf(0.8961);
            let b = -4.0 * ((self.rough / self.ID / 3.7065) - (5.0452 / Re) * (a.log10())).log10();
            1.0 / b.powf(2.0)
        }
    }

    fn AnnularModel(&mut self) {
        let area = std::f64::consts::PI / 4.0 * self.ID * self.ID; // pipe inside cross section area [m^2]
        let UGS = self.WG / (self.LoG * area) / 3600.0; // Superficial Vapor velocity [m/s]
        let ULS = self.WL / (self.LoL * area) / 3600.0; // Superficial Liquid velocity [m/s]

        // Assuming Fanning is a function that you have defined elsewhere
        let fSL = self.fanning(self.LoL * ULS * self.ID / self.muL); // Fanning friction factor for Liquid Phase only in pipe [-]
        let fSG = self.fanning(self.LoG * UGS * self.ID / self.muG); // Fanning friction factor for Vapor Phase only in pipe [-]

        let X2 = fSL * self.LoL * ULS.powi(2) / (fSG * self.LoG * UGS.powi(2)); // Martinelli parameter [-]
        let Y = G * (self.LoL - self.LoG) / (4.0 * fSG * self.LoG * UGS.powi(2) / (2.0 * self.ID)); // Martinelli parameter [-]

        self.alfaL = 0.5; // initial value for Liquid Hold-Up [-]
        let mut delta: f64; // absolute error [-]
        let eps = 1e-4; // allowable tolerance [-]
        let mut gx; // Liquid Holdup function eq.(29)
        let mut gpx; // 1st order derivated function
        let mut alfaLcal; // alfaL (cal.) [-]

        loop {
            gx = X2 * (1.0 - self.alfaL).powf(2.5)
                - self.alfaL.powf(2.0)
                - 75.0 * self.alfaL.powf(3.0)
                - Y * (1.0 - self.alfaL).powf(2.5) * self.alfaL.powf(3.0);

            gpx = -2.5 * X2 * (1.0 - self.alfaL).powf(1.5)
                - 2.0 * self.alfaL
                - 225.0 * self.alfaL.powf(2.0)
                - 3.0 * Y * (1.0 - self.alfaL).powf(2.5) * self.alfaL.powf(2.0)
                + 2.5 * Y * (1.0 - self.alfaL).powf(1.5) * self.alfaL.powf(3.0);

            alfaLcal = self.alfaL - gx / gpx;
            let delta = (self.alfaL - alfaLcal).abs();
            self.alfaL = alfaLcal;

            if delta <= eps {
                break;
            }
        }
        self.Pfric = 2.0 * fSG * self.LoG * UGS.powi(2) / (G * self.ID) * (1.0 + 75.0 * self.alfaL)
            / (1.0 - self.alfaL).powf(2.5)
            / 10000.0
            * 100.0
            * self.SF;
        self.Pgrav = self.LoG / 10000.0 * 100.0;
        self.LoTP = self.LoL * self.alfaL + self.LoG * (1.0 - self.alfaL);
        self.UTP = UGS + ULS;
        let LoNS = (self.WL + self.WG) / (self.WL / self.LoL + self.WG / self.LoG);
        self.Head = LoNS * self.UTP.powi(2) / (2.0 * G) / 10000.0;
        self.Ef = (LoNS * 0.062428) * (self.UTP * 3.28084).powi(2) / 10000.0; // must transfer to imperial unit
    }

    fn SlugModel(&mut self) {
        let area = std::f64::consts::PI / 4.0 * self.ID * self.ID; // pipe inside cross section area [m^2]
        let UGS = self.WG / (self.LoG * area) / 3600.0; // Superficial Vapor velocity [m/s]
        let ULS = self.WL / (self.LoL * area) / 3600.0; // Superficial Liquid velocity [m/s]
        let UTP = UGS + ULS; // Two Phase Velocity [m/s]
        let Um = UTP;
        let C0 = 1.0; // The Distribution parameter [-]
        let K = -0.6; // Drift-flux coefficient [-]
        let Ub = C0 * Um + K * ((G * self.ID * (self.LoL - self.LoG) / self.LoL).sqrt()); // bubble velocity down flow [m/s]
        self.HL = 1.0 - UGS / Ub; // Liquid Hold-up [-]
        if self.HL > 0.75 {
            self.HL = 0.75;
        }
        let alfa = 0.25; // Gas average void fraction [-]
        self.LoLS = self.LoL * (1.0 - alfa) + self.LoG * alfa;
        let muLS = self.muL * (1.0 - alfa) + self.muG * alfa;
        let ReLS = self.LoLS * Um * self.ID / muLS;
        let f0 = self.fanning(ReLS) * 4.0; // single phase Moddy Darcy Friction Factor [-]
        let Landa: f64 = 0.75;
        let LnLanda = -1.0 * Landa.ln();
        let fTP = (1.0
            + LnLanda
                / (1.281 - 0.478 * LnLanda + 0.444 * LnLanda.powf(2.0)
                    - 0.094 * LnLanda.powf(3.0)
                    + 0.00843 * LnLanda.powf(4.0)))
            * f0; // Two Phase Moddy (Darcy) friction factor [-]
        self.Pfric = fTP * self.LoLS * Um.powf(2.0) / (2.0 * G * self.ID) * self.HL / 10000.0
            * 100.0
            * self.SF;
        self.Pgrav = (self.HL * self.LoL + (1.0 - self.HL) * self.LoG) / 10000.0 * 100.0;
        let LoNS = (self.WL + self.WG) / (self.WL / self.LoL + self.WG / self.LoG); // No-Slip Velocity [m/s]
        self.Loip = self.LoL * self.HL + self.LoG * (1.0 - self.HL);
        self.Head = LoNS * UTP.powf(2.0) / (2.0 * G) / 10000.0; // 1.0 Velocity Head
        self.Ef = (LoNS * 0.062428) * (UTP * 3.28084).powf(2.0) / 10000.0; // Erosion Factor must transfer to imperial unit
    }

    fn BubbleModel(&mut self) {
        let area = std::f64::consts::PI / 4.0 * self.ID * self.ID; // pipe inside cross section area [m^2]
        let UGS = self.WG / (self.LoG * area) / 3600.0; // Superficial Vapor velocity [m/s]
        let ULS = self.WL / (self.LoL * area) / 3600.0; // Superficial Liquid velocity [m/s]
        self.UTP = UGS + ULS; // Two Phase Velocity [m/s]
        let Um = self.UTP;
        let C0 = 1.0; // The Distribution parameter [-]
        let K = 0.0; // Drift-flux coefficient [-]
        let Ub = C0 * Um + K * (self.ST * G * (self.LoL - self.LoG) / self.LoL.powi(2)).powf(0.25); // bubble velocity down flow [m/s]
        self.HL = 1.0 - UGS / Ub; // Liquid Hold-up [-]
        let Landa = (self.WL / self.LoL) / (self.WL / self.LoL + self.WG / self.LoG);
        self.LoTP =
            self.LoL * Landa.powi(2) / self.HL + self.LoG * (1.0 - Landa).powi(2) / (1.0 - self.HL); // Two phase density [kg/m^3]
        let muTP = self.muL * Landa + self.muG * (1.0 - Landa); // Two Phase Viscosity [Kg/(m-s)]
        let ReTP = self.LoTP * Um * self.ID / muTP; // Two phase Reynold Number [-]

        // Assuming Fanning is a function that you have defined elsewhere
        let f0 = self.fanning(ReTP) * 4.0; // Darcy friction factor [-]
        let LnLanda = -1.0 * Landa.ln();
        let fTP = (1.0
            + LnLanda
                / (1.281 - 0.478 * LnLanda + 0.444 * LnLanda.powi(2) - 0.094 * LnLanda.powi(3)
                    + 0.00843 * LnLanda.powi(4)))
            * f0; // Two Phase Moddy (Darcy) friction factor [-]
        self.Pfric = fTP * self.LoTP * Um.powi(2) / (2.0 * G * self.ID) / 10000.0 * 100.0 * self.SF;
        self.Pgrav = (self.HL * self.LoL + (1.0 - self.HL) * self.LoG) / 10000.0 * 100.0;
        let LoNS = (self.WL + self.WG) / (self.WL / self.LoL + self.WG / self.LoG); // No-Slip Velocity [m/s]
        self.Head = LoNS * self.UTP.powi(2) / (2.0 * G) / 10000.0; // 1.0 Velocity Head
        self.Ef = (LoNS * 0.062428) * (self.UTP * 3.28084).powi(2) / 10000.0; // Erosion Factor must transfer to imperial unit
    }
}

impl TwoPhaseLine for VerticalDown {
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
        let ratio_a = 0.0;
        let ratio_b = 0.0;
        let ratio_c = 0.0;
        let ratio_d = 0.0;
        let mut Dcrit = 0.0;

        let Vg = self.muG / self.LoG; // vapor kinematic viscosity [m^2/s]
        let VL = self.muL / self.LoL; // Liquid kinematic viscosity [m^2/s]
        let area = std::f64::consts::PI / 4.0 * self.ID * self.ID; // pipe inside cross section area [m^2]
        let UX = self.WG / (self.LoG * area) / 3600.0; // Vapor velocity [m/s]
        let UY = self.WL / (self.LoL * area) / 3600.0; // Liquid velocity [m/s]

        let Db = 0.096887; // 無因次液膜厚度 [-] Db = δ/D
        let SL = std::f64::consts::PI * self.ID; // 氣泡施予管壁之濕潤周長 [m]
        let Si = std::f64::consts::PI * self.ID * (1.0 - 2.0 * Db); // 界面剪應力施予氣液界面的濕潤周長 [m]
        let AL = std::f64::consts::PI * self.ID.powi(2) * (Db - Db.powi(2)); // 管中液膜所佔橫截面積 [m^2]
        let AG = std::f64::consts::PI * self.ID.powi(2) * (0.5 - Db).powi(2); // 管中氣體核所佔橫截面積 [m^2]
        let DL = 4.0 * self.ID * (Db - Db.powi(2)); // 液體的水力直徑 [m]
        let DG = (1.0 - 2.0 * Db) * self.ID; // 氣泡的水力直徑 [m]

        const n: f64 = 0.2; // Eq (7) 中的次冪
        const m: f64 = 0.2; // Eq (7) 中的次冪
        const CG: f64 = 0.046; // 氣體摩擦因子關聯式中的常數
        const CL: f64 = 0.046; // 液體摩擦因子關聯式中的常數
        let UG = 4.0 * UX / (1.0 - 4.0 * Db + 4.0 * Db * Db); // 氣體的真實速度 [m/s]
        let f1 = CL * (DL / VL).powf(-n);
        let f2 = CG * (DG / Vg).powf(-m);
        let fi = f2 * UG.powf(-m);
        let K1 = Si * (1.0 / AL + 1.0 / AG);
        let K2 = fi * self.LoG / 2.0;
        let K3 = K2 * 2.0 * UG;
        let K4 = G * (self.LoL - self.LoG);
        let K5 = f1 * self.LoL / 2.0 * SL / AL;

        // Root-finding by Newton-Raphson Method to solve UL
        let mut delta = 1.0; // absolute error
        let eps = 0.001; // allowable tolerance
        let mut UL = UY / 4.0 / (Db - Db * Db); // initial value of UL

        while delta > eps {
            let gx =
                (K2 * UL.powi(2) - K3 * UL + K2 * UG.powi(2)) * K1 + K4 - K5 * UL.powf(2.0 - n);
            let gpx = (2.0 * K2 * UL - K3) * K1 - K5 * (2.0 - n) * UL.powf(1.0 - n);
            let ULcal = UL - gx / gpx;
            delta = (ULcal - UL).abs();
            UL = (UL + ULcal) / 2.0;
        }

        // ratio A calculation
        let UYA = UL * 4.0 * (Db - Db * Db); // The manuscript here has divided by 2, Eq (9) does not
        let ratio_a = UY / UYA;

        // Dcrit calculation
        Dcrit = 4.36f64.powi(2) * ((self.LoL - self.LoG) * self.ST / self.LoL.powi(2) / G).sqrt();

        // ratio D calculation
        let mut alfa = 0.52;
        let mut U0 = 1.53 * (G * (self.LoL - self.LoG) * self.ST / self.LoL.powi(2)).powf(0.25);
        let UYD = UX * (1.0 - alfa) / alfa + (1.0 - alfa) * U0; // ref. Eq. (19) (24)
        let ratio_d = UY / UYD; // The manuscript here is misprinted as ratioC

        // ratio C calculation
        // Assuming getUYCFromCurveC is a function that you have defined elsewhere
        let UYC = self.get_uyc_from_curve_c(UX);
        let ratio_c = UY / UYC;

        // ratio B calculation
        alfa = 0.25;
        U0 = 1.53 * (G * (self.LoL - self.LoG) * self.ST / self.LoL.powi(2)).powf(0.25);
        let UYB = (1.0 - alfa) / alfa * UX + (1.0 - alfa) * U0; // Eq. (23)
        let ratio_b = UY / UYB;

        // println!("ratio a: {}", ratio_a);
        // println!("ratio b: {}", ratio_b);
        // println!("ratio c: {}", ratio_c);
        // println!("ratio d: {}", ratio_d);
        if ratio_a < 1.0 {
            self.regime_enum = Regime::VerticalDownAnnularFlow(String::from("Annular Flow"));
        } else {
            if self.ID <= Dcrit {
                // Case II, Figure 2(b), Curve C-D
                if ratio_d < 1.0 {
                    self.regime_enum = Regime::VerticalDownSlugFlow(String::from("Slug Flow"));
                } else {
                    if ratio_c < 1.0 {
                        self.regime_enum = Regime::VerticalDownSlugFlow(String::from("Slug Flow"));
                    } else {
                        self.regime_enum = Regime::VerticalDownDispersedBubbleFlow(String::from(
                            "Dispersed-Bubble Flow",
                        ));
                    }
                }
            } else {
                // D > Dcrit , Case I, Figure 2(a), Curve B-C-D
                if ratio_d < 1.0 {
                    self.regime_enum = Regime::VerticalDownSlugFlow(String::from("Slug Flow"));
                } else {
                    if ratio_c < 1.0 {
                        self.regime_enum = Regime::VerticalDownSlugFlow(String::from("Slug Flow"));
                    } else {
                        if ratio_b < 1.0 {
                            self.regime_enum =
                                Regime::VerticalDownSlugFlow(String::from("Slug Flow"));
                        } else {
                            self.regime_enum = Regime::VerticalDownDispersedBubbleFlow(
                                String::from("Dispersed-Bubble Flow"),
                            );
                        }
                    }
                }
            }
        }
        self.flow_regime = match &self.regime_enum {
            Regime::VerticalDownAnnularFlow(v) => v.to_string(),
            Regime::VerticalDownSlugFlow(v) => v.to_string(),
            Regime::VerticalDownDispersedBubbleFlow(v) => v.to_string(),
            _ => String::from(""),
        };
    }

    fn model_cal(&mut self) {
        if !self.is_unit_change {
            self.unit_transfer();
        }
        self.flow_regime();
        match self.regime_enum {
            Regime::VerticalDownAnnularFlow(..) => self.AnnularModel(),
            Regime::VerticalDownSlugFlow(..) => self.SlugModel(),
            Regime::VerticalUpBubbleFlow(..) => self.BubbleModel(),
            Regime::VerticalDownDispersedBubbleFlow(..) => self.BubbleModel(),
            _ => {
                println!("No match model for this flow pattern !!")
            }
        }
    }
}
