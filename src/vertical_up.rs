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
    pub flow_regime: Regime,     // identify the flow regime

    // state variable
    is_unit_change: bool,    // is call the unit_transfer and transfer to unit
}

impl VerticalUp {
    pub fn new(&mut self, WL: f64, WG: f64, LoL: f64, LoG: f64, muL: f64, muG: f64, ST: f64,
               rough: f64, SF: f64, ID: f64, degree: f64) {
        self.WL = WL;
        self.WG = WG;
        self.LoL = LoL;
        self.LoG = LoG;
        self.muL = muL;
        self.muG = muG;
        self.ST = ST;
        self.rough = rough;
        self.SF = SF;
        self.ID = ID;
        self.degree = degree;
        self.is_unit_change = false;
        self.flow_regime = Regime::NONE;
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

    fn flow_regime(&self) -> Regime {
        todo!()
    }

    fn model_cal() {
        todo!()
    }
}