#![allow(dead_code)]
#![allow(unused_variables)]
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]
#![allow(unused_assignments)]

use std::f32::consts::PI;

use twoline::TwoPhaseLine;

use crate::twoline;
use crate::twoline::Regime;

const G: f64 = 9.81;         // gravity accelerator [,/s^2]
const GC: f64 = 9.8;        // gravity constant [kg-m/kgf-s^2]

pub struct VerticalUp {
    // process data
    pub WL: f64,        // liquid mass flow rate [kg/hr]
    pub WG: f64,        // Vapor mass flow rate [kg/hr]
    pub LoL: f64,       // Liquid density [kg/m^3]
    pub LoG: f64,       // Vapor density [kg/m^3]
    pub muL: f64,       // Liquid viscosity [cP] -> [kg/m-s]
    pub muG: f64,       // Vapor viscosity [cP] -> [kg/m-s]
    pub ST: f64,        // Liquid surface tension [dyne/cm] -> [kg/s^2]
    pub rough: f64,     // pipe absolute roughness [mm] -> [m]
    pub SF: f64,        // Safety factor [-]
    pub ID: f64,        // pipe inside diameter [in] -> [m]
    pub degree: f64,    // degree,  Horizontal = 0, -Up / +Down

    // result
    pub regime_enum: Regime,     // identify the flow regime(enum)
    pub flow_regime: String,     // identify the flow regime(String)

    // state variable
    is_unit_change: bool,    // is call the unit_transfer and transfer to unit
}

impl VerticalUp {
    pub fn new(WL: f64, WG: f64, LoL: f64, LoG: f64, muL: f64, muG: f64, ST: f64,
               rough: f64, SF: f64, ID: f64, degree: f64) -> Self {
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
        }
    }

    fn get_UGSE_from_curveE(&self) -> f64 {
        // Refer to Eq. (21)
        let term_e = ((self.LoL - self.LoG) * G * self.ST).powf(0.25);
        let UGS_cal = 3.1 * term_e / self.LoG.sqrt();      // Curve E 求得的 UGS 計算值 (Curve E 與 ULS 無關)
        UGS_cal
    }

    fn get_ULSB_from_curveB(&self, x: f64) -> f64 {
        // by Eq. (12)
        let term_b = 4.0 * ((G * (self.LoL - self.LoG) / self.LoL).powf(0.446) * self.ID.powf(0.429)
            * (self.ST / self.LoL).powf(0.089) / (self.muL / self.LoL).powf(0.072));
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
}

impl TwoPhaseLine for VerticalUp {
    fn unit_transfer(&mut self) {
        if !self.is_unit_change {
            self.muL = self.muL * 0.001;    // [cP] -> [kg/m-s]
            self.muG = self.muG * 0.001;    // [cP] -> [kg/m-s]
            self.ID = self.ID * 2.54 / 100.0;       // [in] -> [m]
            self.ST = self.ST * 1.019716213E-4;     // [dyne/cm] -> [kgf/s^2]
            self.degree = self.degree * (PI as f64) / 180.0;        // [degree] -> [rad]
            self.rough = self.rough / 1000.0;       // [mm] -> [m]

            self.is_unit_change = true;
        }
    }

    fn flow_regime(&mut self) {
        let alfa = 0.25;                     // Average Gas Void Fraction
        let area = std::f64::consts::PI * self.ID * self.ID / 4.0;    // pipe area [m^2]
        let UG = self.WG / self.LoG / area / 3600.0;                  // Vapor Velocity [m/s]
        let UL = self.WL / self.LoL / area / 3600.0;                  // Liquid Velocity [m/s]
        let UGS = UG * alfa;                 // Vapor Superficial Velocity [m/s]
        let ULS = UL * (1.0 - alfa);         // Liquid Superficial Velocity [m/s]

        // Curve E
        let ratio_e = UGS / self.get_UGSE_from_curveE();

        // Curve C
        let ratio_c = UGS / (13.0 / 12.0 * ULS);    // Curve C Eq. (15) 求得的 UGS 計算值

        // Curve B
        let ratio_b = ULS / self.get_ULSB_from_curveB(UGS);

        // Curve A ()
        let mut ratio_a: f64 = 0.0;
        match self.get_UGSA_from_curveA(ULS) {
            Ok(value) => { ratio_a = UGS / value; }
            Err(e) => { println!("Error: {}", e); }
        }

        // ***** Regime 的判斷邏輯 *****
        if ratio_e > 1.0 {
            // Churn transition to Annular Flow 與流體速度無關, 與管徑亦無任何關聯
            // ratioE > 1 : Annular Flow
            // ratioE <= 1 : Churn Flow
            self.regime_enum = Regime::VerticalUpAnnularFlow(String::from("Vertical Up Annular Flow"));
        } else if ratio_a <= 1.0 && ratio_b <= 1.0 {
            self.regime_enum = Regime::VerticalUpBubbleFlow(String::from("Vertical Up Bubble Flow"));
        } else if ratio_a > 1.0 && ratio_b <= 1.0 {
            self.regime_enum = Regime::VerticalUpSlugAndChurnFlow(String::from("Vertical Up Slug and Churn Flow"));
        } else if ratio_a > 1.0 && ratio_c > 1.0 {
            self.regime_enum = Regime::VerticalUpSlugAndChurnFlow(String::from("Vertical Up Slug and Churn Flow"));
        } else {
            self.regime_enum = Regime::VerticalUpFinelyDispersedBubbleFlow(String::from("Vertical Up Finely Dispersed Bubble Flow"));
        };

        self.flow_regime = match &self.regime_enum {
            Regime::VerticalUpAnnularFlow(v) => v.to_string(),
            Regime::VerticalUpBubbleFlow(v) => v.to_string(),
            Regime::VerticalUpSlugAndChurnFlow(v) => v.to_string(),
            Regime::VerticalUpFinelyDispersedBubbleFlow(v) => v.to_string(),
            _ => String::from(""),
        };
    }

    #[allow(dead_code)]
    fn model_cal() {
        todo!()
    }
}