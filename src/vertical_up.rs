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

pub struct VerticalUp {
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

    // for Similarity Analysis Model
    pub Loip: f64,  // two phase density [kg/m^3]
    pub RL: f64,    // Liquid Volume Fraction [-]
    pub UTP: f64,   // Two-Phase Velocity [m/sec]
    pub Head: f64,  // 1.0 Velocity Head [kgf/cm^2]
    pub Pfric: f64, // Frictional Pressure Loss [kgf/cm^2/100m]
    pub Pgrav: f64, // Elevation Head Loss [kgf/cm^2/100m]
    pub Ef: f64,    // Errosion Factor [-]

    // for Bubble Model
    pub LoNS: f64,  // two phase density [kg/m^3]
    pub Landa: f64, // Liquid Volume Fraction [-]
    // UTP, Head, Pfric, Pgrav, Ef same as Similarity

    // for Slug Model
    pub LoLS: f64, // Liquid Slug Density [kg/m^3]
    pub LoSU: f64, // Two phase slug unit density [kg/m^3]
    pub ULLS: f64, // Liquid Slug Velocity [m/sec]
    pub LLS: f64,  // Liquid Slug Length [m]
    pub Lu: f64,   // Slug unit length [m]
    pub Le: f64,   // Stabilizes to Slug flow in x m
    // Head, Pfric, Pgrav, Ef same as Similarity

    // state variable
    is_unit_change: bool, // is call the unit_transfer and transfer to unit
}

impl VerticalUp {
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
        VerticalUp {
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
            Loip: 0.0,
            RL: 0.0,
            UTP: 0.0,
            Head: 0.0,
            Pfric: 0.0,
            Pgrav: 0.0,
            Ef: 0.0,
            LoNS: 0.0,
            Landa: 0.0,
            LoLS: 0.0,
            LoSU: 0.0,
            ULLS: 0.0,
            LLS: 0.0,
            Lu: 0.0,
            Le: 0.0,
        }
    }

    fn get_UGSE_from_curveE(&self) -> f64 {
        // Refer to Eq. (21)
        let term_e = ((self.LoL - self.LoG) * G * self.ST).powf(0.25);
        let UGS_cal = 3.1 * term_e / self.LoG.sqrt(); // Curve E 求得的 UGS 計算值 (Curve E 與 ULS 無關)
        UGS_cal
    }

    fn get_ULSB_from_curveB(&self, x: f64) -> f64 {
        // by Eq. (12)
        let term_b = 4.0
            * ((G * (self.LoL - self.LoG) / self.LoL).powf(0.446)
                * self.ID.powf(0.429)
                * (self.ST / self.LoL).powf(0.089)
                / (self.muL / self.LoL).powf(0.072));
        let ULS_cal = term_b - x;
        ULS_cal
    }

    fn get_UGSA_from_curveA(&self, y: f64) -> Result<f64, &'static str> {
        // by Eq. (5)
        let term_a = G * (self.LoL - self.LoG) * self.ST / self.LoL.powf(2.0);
        let UGS_cal = (y + 0.9938 * term_a.powf(0.25)) / 3.0;

        if UGS_cal.is_nan() || UGS_cal.is_infinite() {
            Err("VerticalUp-DT Curve A: ArithmeticException")
        } else {
            Ok(UGS_cal)
        }
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

    fn SimilarityAnalysis(&mut self) {
        // for Anaular flow pattern
        use std::f64;
        let area = f64::consts::PI * self.ID * self.ID / 4.0; // pipe area [m^2]
        let Gt = (self.WL + self.WG) / area / 3600.0; // Eq (22)

        let UGS = self.WG / self.LoG / area / 3600.0; // Vapor Velocity [m/s]
        let ULS = self.WL / self.LoL / area / 3600.0; // Liquid Velocity [m/s]
        self.UTP = UGS + ULS; // Two Phase Velocity [m/s], Eq (23)

        let lamda = ULS / (ULS + UGS); // Liquid Volume Fraction [-], Eq (24)
        let mut Rgi = 0.5; // Gas Hold-up (Rg) initial value [-]
        let eps = 1e-4; // allowable tolerance
        let np = 100; // trial number
        let i = 0; // loop counter

        for i in 0..np {
            // (5) Calc. Re and Fr
            let Re = self.ID * Gt / (Rgi * self.muG + (1.0 - Rgi) * self.muL); // Eq. (25)
            let Fr = self.UTP * self.UTP / (G * self.ID); // Froude Number, Eq. (26)

            // (6) Calc. Z and K
            let Z = Re.powf(0.167) * Fr.powf(0.125) / lamda.powf(0.25); // Eq.(27)
            let K;
            if Z < 10.0 {
                K = -0.16367 + 0.31037 * Z - 0.03525 * Z * Z + 0.001366 * Z * Z * Z;
            } else {
                K = 0.75545 + 0.003585 * Z - 0.1436e-4 * Z * Z;
            }

            // (7) Calc Rg (cal.)
            let x = self.WG / (self.WG + self.WL);
            let Rgcal = K / ((1.0 / x - 1.0) * (self.LoG / self.LoL) + 1.0); // Eq. (28)

            // (8) Calc delta and judgement convergence condition
            let delta = (Rgcal - Rgi).abs();
            if delta > eps {
                Rgi = (Rgcal + Rgi) / 2.0;
                // Repeat calc (5), (6), (7)
            } else {
                break;
            }
        }

        let Rg;
        if i < np {
            Rg = Rgi; // certain Rg
        } else {
            Rg = 0.0; // no convergence
            return;
        }

        // Calculate Result
        self.RL = 1.0 - Rg;
        let LoTP =
            self.LoL * lamda.powf(2.0) / (1.0 - Rg) + self.LoG * (1.0 - lamda).powf(2.0) / Rg; // Eq. (29)
        let muTP = self.muL * lamda + self.muG * (1.0 - lamda);
        let ReTP = self.ID * (ULS + UGS) * LoTP / muTP; // Eq. (30)
        let f0 = self.fanning(ReTP) * 4.0;
        let LnLanda = -1.0 * lamda.ln();
        let fTP = (1.0
            + LnLanda
                / (1.281 - 0.478 * LnLanda + 0.444 * LnLanda.powf(2.0)
                    - 0.094 * LnLanda.powf(3.0)
                    + 0.00843 * LnLanda.powf(4.0)))
            * f0; // Eq. (31)

        self.Pfric =
            fTP * LoTP * (ULS + UGS).powf(2.0) / (2.0 * G * self.ID) / 10000.0 * 100.0 * self.SF; // Eq. (32)
        self.Pgrav = (self.LoL * (1.0 - Rg) + self.LoG * Rg) / 10000.0 * 100.0; // Eq. (33)
        self.Loip = self.LoL * (1.0 - Rg) + self.LoG * Rg;
        let LoNS = (self.WL + self.WG) / (self.WL / self.LoL + self.WG / self.LoG);
        self.Head = LoNS * self.UTP.powf(2.0) / (2.0 * G) / 10000.0;
        self.Ef = (LoNS * 0.062428) * ((ULS + UGS) * 3.28084).powf(2.0) / 10000.0;
    }

    fn SlugModel(&mut self) {
        // for Slug and Churn flow pattern
        use std::f64;

        let area = f64::consts::PI * self.ID * self.ID / 4.0; // pipe area [m^2]
        let UGS = self.WG / self.LoG / area / 3600.0; // Vapor Velocity [m/s]
        let ULS = self.WL / self.LoL / area / 3600.0; // Liquid Velocity [m/s]
        let UTP = UGS + ULS; // Two Phase Velocity [m/s]

        let alfaLS = 0.25; // Void fraction of liquid slug
        let UN = 0.35 * (G * self.ID).sqrt() + 1.29 * UTP; // rise velocity of transition of Taylor Bubble in stagnant liquid [m/s], Eq. (A-16)

        // term = UO : the rise velocity due to buoyancy [m/s], Eq. (A-18)
        let mut term;
        term = 1.53
            * ((self.ST * G * (self.LoL - self.LoG)) / (self.LoL * self.LoL)).powf(0.25)
            * (1.0f64 - alfaLS).sqrt();
        self.ULLS = UTP - term * alfaLS; // Velocity of the liquid in the liquid slug
        let UGLS = UTP + term * (1.0 - alfaLS); // Velocity of the gas in the liquid slug
        let Landa = self.ULLS * (1.0 - alfaLS) / UTP; // Liquid volume fraction [-]
        let mut alfaTB = 0.1; // Void fraction of Taylor Bubble [-]
        let mut delta = 1.0; // absolute error
        let eps = 1e-4; // allowable tolerance

        let mut term1;
        let mut alf1;
        let mut alf2;
        let mut alfaTB_cal;

        while delta > eps {
            term = UN * (alfaTB - alfaLS) - self.ULLS * (1.0 - alfaLS);
            term1 = 9.916 * (G * self.ID * (1.0 - alfaTB.sqrt())).sqrt() * (1.0 - alfaTB);
            alf1 = term - term1;
            alf2 = UN
                + 9.916
                    * (G * self.ID * (1.0 - alfaTB.sqrt())).sqrt()
                    * (1.0 + 5.0 * alfaTB.sqrt())
                    / (4.0 * alfaTB.sqrt());
            if alf2 == 0.0 {
                break;
            }
            alfaTB_cal = alfaTB - alf1 / alf2; // Calc. Void fraction of Taylor Bubble [-]
            delta = (alfaTB_cal - alfaTB).abs();
            alfaTB = (alfaTB + alfaTB_cal) / 2.0;
        }

        let ULTB = 9.916 * (G * self.ID * (1.0 - alfaTB.sqrt())).sqrt();
        let beta = (UGS - alfaLS * UGLS) / UN / (alfaTB - alfaLS); // LTB/Lu
        let alfaSU = beta * alfaTB + (1.0 - beta) * alfaLS; // void fraction of a slug unit
        self.LLS = 20.0 * self.ID;
        self.Lu = self.LLS / (1.0 - beta);
        let LTB = self.Lu - self.LLS; // length of Taylor Bubble
        self.LoLS = self.LoG * (1.0 - Landa).powf(2.0) / alfaLS
            + self.LoL * Landa.powf(2.0) / (1.0 - alfaLS);
        self.LoSU = self.LoG * (1.0 - alfaSU) + self.LoL * alfaSU;
        self.Le = self.ID * 35.5 * (8.0 / 7.0 * UTP / (G * self.ID).sqrt() + 0.25) * 1.2;

        let LoNS = (self.WL + self.WG) / (self.WL / self.LoL + self.WG / self.LoG); // no-slip density [Kg/m^3]
        self.Head = LoNS * UTP.powf(2.0) / (2.0 * G) / 10000.0;
        let LoTP = self.LoLS;
        let muTP = self.muL * Landa + self.muG * (1.0 - Landa);
        let ReTP = LoTP * UTP * self.ID / muTP;
        let f0 = self.fanning(ReTP) * 4.0;
        let LnLanda = -1.0 * Landa.ln();
        let fTP = f0
            * (1.0
                + LnLanda
                    / (1.281 - 0.478 * LnLanda + 0.444 * LnLanda.powf(2.0)
                        - 0.094 * LnLanda.powf(3.0)
                        + 0.00843 * LnLanda.powf(4.0)));

        self.Pfric = fTP * LoTP * self.ULLS.powf(2.0) / (2.0 * G * self.ID) * (self.LLS / self.Lu)
            / 10000.0
            * 100.0
            * self.SF;
        let Pacc = self.LoL * ULTB / G * (1.0 - alfaTB) * (self.ULLS + ULTB) * (1.0 / self.Lu)
            / 10000.0
            * 100.0;
        self.Pfric = self.Pfric + Pacc;
        self.Pgrav = (self.LoL * (1.0 - alfaSU) + self.LoG * alfaSU) / 10000.0 * 100.0;
        self.Ef = (LoNS * 0.062428) * ((ULS + UGS) * 3.28084).powf(2.0) / 10000.0;
        // must transfer to imperial unit
    }

    fn BubbleModel(&mut self) {
        // for Bubble flow and Finely Bubble flow pattern
        use std::f64;

        let area = f64::consts::PI * self.ID * self.ID / 4.0; // pipe area [m^2]
        let UGS = self.WG / self.LoG / area / 3600.0; // Vapor Velocity [m/s]
        let ULS = self.WL / self.LoL / area / 3600.0; // Liquid Velocity [m/s]
        self.Landa = ULS / (ULS + UGS); // Liquid Volume Fraction [-]
        let mut delta = 1.0; // Absolute error [-]
        let mut alfa = 0.5; // Gas average void fraction [-]
        let eps = 1e-4; // Allowable Tolerance
        let mut n = 0; // trials
        const MAX_TRIAL: i32 = 100;

        while delta > eps {
            let U0 = (1.0f64 - alfa).sqrt()
                * 1.53
                * ((self.ST * G * (self.LoL - self.LoG)) / (self.LoL * self.LoL)).powf(0.25); // Eq. (38)
            let alfacal = UGS / (ULS / (1.0 - alfa) + U0); // Eq. (37)
            delta = (alfacal - alfa).abs();
            alfa = (alfa + alfacal) / 2.0;
            n += 1;
            if n > MAX_TRIAL {
                break;
            }
        }

        if n > MAX_TRIAL {
            println!("Vertical-Up Bubble Model: alfa did not converge!");
        }

        // calculate result here
        let loTP = self.LoG * (1.0 - self.Landa).powf(2.0) / alfa
            + self.LoL * self.Landa.powf(2.0) / (1.0 - alfa); // Eq. (40)
        let muTP = self.muL * self.Landa + self.muG * (1.0 - self.Landa); // Eq. (40)
        let ReTP = self.ID * (ULS + UGS) * loTP / muTP; // Eq. (41)
        let f0 = self.fanning(ReTP) * 4.0; // Step (3)
        let lnlanda = -1.0 * self.Landa.ln();
        let fTP = f0
            * (1.0
                + lnlanda
                    / (1.281 - 0.478 * lnlanda + 0.444 * lnlanda.powf(2.0)
                        - 0.094 * lnlanda.powf(3.0)
                        + 0.00843 * lnlanda.powf(4.0))); // Step (4)
        self.Pfric =
            fTP * loTP * (ULS + UGS).powf(2.0) / (2.0 * G * self.ID) / 10000.0 * 100.0 * self.SF; // Eq. (42)
        self.Pgrav = (self.LoL * (1.0 - alfa) + self.LoG * alfa) / 10000.0 * 100.0; // Eq. (43)
        self.UTP = ULS + UGS;
        self.LoNS = self.LoL * self.Landa + self.LoG * (1.0 - self.Landa);
        self.Head = self.LoNS * self.UTP.powf(2.0) / (2.0 * G) / 10000.0;
        self.Ef = (self.LoNS * 0.062428) * (self.UTP * 3.28084).powf(2.0) / 10000.0;
        // must transfer to imperial unit
    }
}

impl TwoPhaseLine for VerticalUp {
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
        let alfa = 0.25; // Average Gas Void Fraction
        let area = std::f64::consts::PI * self.ID * self.ID / 4.0; // pipe area [m^2]
        let UG = self.WG / self.LoG / area / 3600.0; // Vapor Velocity [m/s]
        let UL = self.WL / self.LoL / area / 3600.0; // Liquid Velocity [m/s]
        let UGS = UG * alfa; // Vapor Superficial Velocity [m/s]
        let ULS = UL * (1.0 - alfa); // Liquid Superficial Velocity [m/s]

        // Curve E
        let ratio_e = UGS / self.get_UGSE_from_curveE();

        // Curve C
        let ratio_c = UGS / (13.0 / 12.0 * ULS); // Curve C Eq. (15) 求得的 UGS 計算值

        // Curve B
        let ratio_b = ULS / self.get_ULSB_from_curveB(UGS);

        // Curve A ()
        let mut ratio_a: f64 = 0.0;
        match self.get_UGSA_from_curveA(ULS) {
            Ok(value) => {
                ratio_a = UGS / value;
            }
            Err(e) => {
                println!("Error: {}", e);
            }
        }

        // ***** Regime 的判斷邏輯 *****
        if ratio_e > 1.0 {
            // Churn transition to Annular Flow 與流體速度無關, 與管徑亦無任何關聯
            // ratioE > 1 : Annular Flow
            // ratioE <= 1 : Churn Flow
            self.regime_enum =
                Regime::VerticalUpAnnularFlow(String::from("Vertical Up Annular Flow"));
        } else if ratio_a <= 1.0 && ratio_b <= 1.0 {
            self.regime_enum =
                Regime::VerticalUpBubbleFlow(String::from("Vertical Up Bubble Flow"));
        } else if ratio_a > 1.0 && ratio_b <= 1.0 {
            self.regime_enum =
                Regime::VerticalUpSlugAndChurnFlow(String::from("Vertical Up Slug and Churn Flow"));
        } else if ratio_a > 1.0 && ratio_c > 1.0 {
            self.regime_enum =
                Regime::VerticalUpSlugAndChurnFlow(String::from("Vertical Up Slug and Churn Flow"));
        } else {
            self.regime_enum = Regime::VerticalUpFinelyDispersedBubbleFlow(String::from(
                "Vertical Up Finely Dispersed Bubble Flow",
            ));
        };

        self.flow_regime = match &self.regime_enum {
            Regime::VerticalUpAnnularFlow(v) => v.to_string(),
            Regime::VerticalUpBubbleFlow(v) => v.to_string(),
            Regime::VerticalUpSlugAndChurnFlow(v) => v.to_string(),
            Regime::VerticalUpFinelyDispersedBubbleFlow(v) => v.to_string(),
            _ => String::from(""),
        };
    }

    fn model_cal(&mut self) {
        if !self.is_unit_change {
            self.unit_transfer();
        }
        self.flow_regime();
        match self.regime_enum {
            Regime::VerticalUpAnnularFlow(..) => self.SimilarityAnalysis(),
            Regime::VerticalUpSlugAndChurnFlow(..) => self.SlugModel(),
            Regime::VerticalUpBubbleFlow(..) => self.BubbleModel(),
            Regime::VerticalUpFinelyDispersedBubbleFlow(..) => self.BubbleModel(),
            _ => {
                println!("No match model for this flow pattern !!")
            }
        }
    }
}
