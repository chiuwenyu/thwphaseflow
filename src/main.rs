use crate::twoline::TwoPhaseLine;

mod vertical_up;
mod twoline;
mod horizontal;

fn main() {
    //Region Test data for Annular Flow
    // // Liquid data
    // let wl: f64 = 72036.365;            // [kg/hr]
    // let lo_l: f64 = 379.63758;           // [kg/m^3]
    // let mu_l: f64 = 0.054;               // [cP]
    // let surface_tension: f64 = 40.0;   // [dyne/cm]
    // // Vapor data
    // let wg: f64 = 78722.747;           // [kg/hr]
    // let lo_g: f64 = 75.286778;           // [kg/m^3]
    // let mu_g: f64 = 0.011;               // [cP]
    // // Misc. data
    // let id: f64 = 12.0;                 // [in]
    // let slope: f64 = 0.0;               // [degree]
    // let rough: f64 = 0.04572;           // [mm]
    // let sf: f64 = 1.0;                  // [-]
    //
    // let mut p1 = VerticalUp::new(wl, wg, lo_l, lo_g, mu_l, mu_g, surface_tension, rough, sf, id, slope);
    // p1.unit_transfer();
    // p1.flow_regime();
    // println!("flow regime << {} >>", p1.flow_regime);
    // p1.model_cal();
    // println!("Tow-Phase Density (kg/m^3) = {:.4}", p1.Loip);
    // println!("Liquid Volume Fraction (-) = {:.3}", p1.RL);
    // println!("Two-Phase Velocity (m/sec) = {:.4}", p1.UTP);
    // println!("1.0 Velocity Head (kgf/cm^2) = {:.4}", p1.Head);
    // println!("Frictional Pressure Loss (kgf/cm^2/100m) = {:.4}", p1.Pfric);
    // println!("Elevation Head Loss (kgf/cm^2/100m) = {:.4}", p1.Pgrav);
    // println!("Erosion Factor (-) = {:.3}", p1.Ef);
    // println!("if Φ ≤ 1 : No Erosion, Φ > 1, Erosion occurred");
    //EndRegion

    // Region Test data for Bubble Flow
    // Liquid data
    // let wl: f64 = 100000.0;            // [kg/hr]
    // let lo_l: f64 = 500.0;           // [kg/m^3]
    // let mu_l: f64 = 1.0;               // [cP]
    // let surface_tension: f64 = 30.0;   // [dyne/cm]
    // // Vapor data
    // let wg: f64 = 50.0;           // [kg/hr]
    // let lo_g: f64 = 2.0;           // [kg/m^3]
    // let mu_g: f64 = 0.01;               // [cP]
    // // Misc. data
    // let id: f64 = 6.065;                 // [in]
    // let slope: f64 = 0.0;               // [degree]
    // let rough: f64 = 0.046;           // [mm]
    // let sf: f64 = 1.0;                  // [-]
    //
    // let mut p2 = VerticalUp::new(wl, wg, lo_l, lo_g, mu_l, mu_g, surface_tension, rough, sf, id, slope);
    // p2.unit_transfer();
    // p2.flow_regime();
    // println!("flow regime << {} >>", p2.flow_regime);
    // p2.model_cal();
    // println!("Tow-Phase Density (kg/m^3) = {:.4}", p2.LoNS);
    // println!("Liquid Volume Fraction (-) = {:.4}", p2.Landa);
    // println!("Two-Phase Velocity (m/sec) = {:.4}", p2.UTP);
    // println!("1.0 Velocity Head (kgf/cm^2) = {:.4}", p2.Head);
    // println!("Frictional Pressure Loss (kgf/cm^2/100m) = {:.4}", p2.Pfric);
    // println!("Elevation Head Loss (kgf/cm^2/100m) = {:.4}", p2.Pgrav);
    // println!("Erosion Factor (-) = {:.3}", p2.Ef);
    // println!("if Φ ≤ 1 : No Erosion, Φ > 1, Erosion occurred");
    //EndRegion

    // Region Test data for Slug Model
    // Liquid data
    // let wl: f64 = 90718.0;            // [kg/hr]
    // let lo_l: f64 = 640.73852;           // [kg/m^3]
    // let mu_l: f64 = 0.3;               // [cP]
    // let surface_tension: f64 = 20.0;   // [dyne/cm]
    // // Vapor data
    // let wg: f64 = 1814.36;           // [kg/hr]
    // let lo_g: f64 = 8.00923;           // [kg/m^3]
    // let mu_g: f64 = 0.01;               // [cP]
    // // Misc. data
    // let id: f64 = 6.065;                 // [in]
    // let slope: f64 = 0.0;               // [degree]
    // let rough: f64 = 0.04572;           // [mm]
    // let sf: f64 = 1.0;                  // [-]
    //
    // let mut p3 = VerticalUp::new(wl, wg, lo_l, lo_g, mu_l, mu_g, surface_tension, rough, sf, id, slope);
    // p3.unit_transfer();
    // p3.flow_regime();
    // println!("flow regime << {} >>", p3.flow_regime);
    // p3.model_cal();
    // println!("Liquid Slug Density (kg/m^3) = {:.4}", p3.LoLS);
    // println!("Tow-Phase Slug Unit Density (kg/m^3) = {:.4}", p3.LoSU);
    // println!("Liquid Slug Velocity (m/sec) = {:.4}", p3.ULLS);
    // println!("Liquid Slug Length (m) = {:.4}", p3.LLS);
    // println!("Slug Unit Length (Liq + Vap) (m) = {:.4}", p3.Lu);
    // println!("Stabilizes to Slug Flow in x m = {:.4}", p3.Le);
    // println!("1.0 Velocity Head (kgf/cm^2) = {:.4}", p3.Head);
    // println!("Frictional Pressure Loss (kgf/cm^2/100m) = {:.4}", p3.Pfric);
    // println!("Elevation Head Loss (kgf/cm^2/100m) = {:.4}", p3.Pgrav);
    // println!("Erosion Factor (-) = {:.3}", p3.Ef);
    // println!("if Φ ≤ 1 : No Erosion, Φ > 1, Erosion occurred");
    //EndRegion
}
