use crate::twoline::TwoPhaseLine;
use crate::vertical_up::VerticalUp;

mod vertical_up;
mod twoline;

fn main() {
    // test data
    // Liquid data
    let wl: f64 = 72036.365;            // [kg/hr]
    let lo_l: f64 = 379.63758;           // [kg/m^3]
    let mu_l: f64 = 0.054;               // [cP]
    let surface_tension: f64 = 40.0;   // [dyne/cm]
    // Vapor data
    let wg: f64 = 78722.747;           // [kg/hr]
    let lo_g: f64 = 75.286778;           // [kg/m^3]
    let mu_g: f64 = 0.011;               // [cP]
    // Misc. data
    let id: f64 = 12.0;                 // [in]
    let slope: f64 = 0.0;               // [degree]
    let rough: f64 = 0.04572;           // [mm]
    let sf: f64 = 1.0;                  // [-]

    let mut p1 = VerticalUp::new(wl, wg, lo_l, lo_g, mu_l, mu_g, surface_tension, rough, sf, id, slope);
    p1.unit_transfer();
    p1.flow_regime();
    println!("flow regime << {} >>", p1.flow_regime);
    p1.model_cal();
    println!("Tow-Phase Density (kg/m^3) = {:.4}", p1.Loip);
    println!("Liquid Volume Fraction (-) = {:.3}", p1.RL);
    println!("Two-Phase Velocity (m/sec) = {:.4}", p1.UTP);
    println!("1.0 Velocity Head (kgf/cm^2) = {:.4}", p1.Head);
    println!("Frictional Pressure Loss (kgf/cm^2/100m) = {:.4}", p1.Pfric);
    println!("Elevation Head Loss (kgf/cm^2/100m) = {:.4}", p1.Pgrav);
    println!("Erosion Factor (-) = {:.3}", p1.Ef);
    println!("if Φ ≤ 1 : No Erosion, Φ > 1, Erosion occurred");
}
