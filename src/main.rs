use crate::twoline::TwoPhaseLine;
use crate::vertical_down::VerticalDown;
use crate::vertical_up::VerticalUp;

mod horizontal;
mod twoline;
mod vertical_down;
mod vertical_up;

pub fn vertical_up_validate() {
    //Region Test data for Annular Flow
    // Liquid data
    let wl: f64 = 72036.365; // [kg/hr]
    let lo_l: f64 = 379.63758; // [kg/m^3]
    let mu_l: f64 = 0.054; // [cP]
    let surface_tension: f64 = 40.0; // [dyne/cm]
                                     // Vapor data
    let wg: f64 = 78722.747; // [kg/hr]
    let lo_g: f64 = 75.286778; // [kg/m^3]
    let mu_g: f64 = 0.011; // [cP]
                           // Misc. data
    let id: f64 = 12.0; // [in]
    let slope: f64 = 0.0; // [degree]
    let rough: f64 = 0.04572; // [mm]
    let sf: f64 = 1.0; // [-]

    let mut p1 = VerticalUp::new(
        wl,
        wg,
        lo_l,
        lo_g,
        mu_l,
        mu_g,
        surface_tension,
        rough,
        sf,
        id,
        slope,
    );
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
    // EndRegion

    // Region Test data for Bubble Flow
    // Liquid data
    let wl: f64 = 100000.0; // [kg/hr]
    let lo_l: f64 = 500.0; // [kg/m^3]
    let mu_l: f64 = 1.0; // [cP]
    let surface_tension: f64 = 30.0; // [dyne/cm]
                                     // Vapor data
    let wg: f64 = 50.0; // [kg/hr]
    let lo_g: f64 = 2.0; // [kg/m^3]
    let mu_g: f64 = 0.01; // [cP]
                          // Misc. data
    let id: f64 = 6.065; // [in]
    let slope: f64 = 0.0; // [degree]
    let rough: f64 = 0.046; // [mm]
    let sf: f64 = 1.0; // [-]

    let mut p2 = VerticalUp::new(
        wl,
        wg,
        lo_l,
        lo_g,
        mu_l,
        mu_g,
        surface_tension,
        rough,
        sf,
        id,
        slope,
    );
    p2.unit_transfer();
    p2.flow_regime();
    println!("flow regime << {} >>", p2.flow_regime);
    p2.model_cal();
    println!("Tow-Phase Density (kg/m^3) = {:.4}", p2.LoNS);
    println!("Liquid Volume Fraction (-) = {:.4}", p2.Landa);
    println!("Two-Phase Velocity (m/sec) = {:.4}", p2.UTP);
    println!("1.0 Velocity Head (kgf/cm^2) = {:.4}", p2.Head);
    println!("Frictional Pressure Loss (kgf/cm^2/100m) = {:.4}", p2.Pfric);
    println!("Elevation Head Loss (kgf/cm^2/100m) = {:.4}", p2.Pgrav);
    println!("Erosion Factor (-) = {:.3}", p2.Ef);
    println!("if Φ ≤ 1 : No Erosion, Φ > 1, Erosion occurred");
    //EndRegion

    // Region Test data for Slug Model
    // Liquid data
    let wl: f64 = 90718.0; // [kg/hr]
    let lo_l: f64 = 640.73852; // [kg/m^3]
    let mu_l: f64 = 0.3; // [cP]
    let surface_tension: f64 = 20.0; // [dyne/cm]
                                     // Vapor data
    let wg: f64 = 1814.36; // [kg/hr]
    let lo_g: f64 = 8.00923; // [kg/m^3]
    let mu_g: f64 = 0.01; // [cP]
                          // Misc. data
    let id: f64 = 6.065; // [in]
    let slope: f64 = 0.0; // [degree]
    let rough: f64 = 0.04572; // [mm]
    let sf: f64 = 1.0; // [-]

    let mut p3 = VerticalUp::new(
        wl,
        wg,
        lo_l,
        lo_g,
        mu_l,
        mu_g,
        surface_tension,
        rough,
        sf,
        id,
        slope,
    );
    p3.unit_transfer();
    p3.flow_regime();
    println!("flow regime << {} >>", p3.flow_regime);
    p3.model_cal();
    println!("Liquid Slug Density (kg/m^3) = {:.4}", p3.LoLS);
    println!("Tow-Phase Slug Unit Density (kg/m^3) = {:.4}", p3.LoSU);
    println!("Liquid Slug Velocity (m/sec) = {:.4}", p3.ULLS);
    println!("Liquid Slug Length (m) = {:.4}", p3.LLS);
    println!("Slug Unit Length (Liq + Vap) (m) = {:.4}", p3.Lu);
    println!("Stabilizes to Slug Flow in x m = {:.4}", p3.Le);
    println!("1.0 Velocity Head (kgf/cm^2) = {:.4}", p3.Head);
    println!("Frictional Pressure Loss (kgf/cm^2/100m) = {:.4}", p3.Pfric);
    println!("Elevation Head Loss (kgf/cm^2/100m) = {:.4}", p3.Pgrav);
    println!("Erosion Factor (-) = {:.3}", p3.Ef);
    println!("if Φ ≤ 1 : No Erosion, Φ > 1, Erosion occurred");
    //EndRegion
}

pub fn horizontal_validate() {
    //Region Test data for Annular-Dispersed Flow (Similarity Model)
    // Liquid data
    // let wl: f64 = 64870.62744; // [kg/hr]
    // let lo_l: f64 = 790.9917; // [kg/m^3]
    // let mu_l: f64 = 0.241; // [cP]
    // let surface_tension: f64 = 14.78; // [dyne/cm]
    //                                   // Vapor data
    // let wg: f64 = 21623.54248; // [kg/hr]
    // let lo_g: f64 = 4.58128; // [kg/m^3]
    // let mu_g: f64 = 0.0091; // [cP]
    //                         // Misc. data
    // let id: f64 = 7.981; // [in]
    // let slope: f64 = 0.0; // [degree]
    // let rough: f64 = 0.04572; // [mm]
    // let sf: f64 = 1.0; // [-]
    //
    // let mut p1 = Horizontal::new(
    //     wl,
    //     wg,
    //     lo_l,
    //     lo_g,
    //     mu_l,
    //     mu_g,
    //     surface_tension,
    //     rough,
    //     sf,
    //     id,
    //     slope,
    // );
    // p1.unit_transfer();
    // p1.flow_regime();
    // println!("p1 flow regime << {} >>", p1.flow_regime);
    // p1.model_cal();
    // println!("Tow-Phase Density (kg/m^3) = {:.4}", p1.Loip);
    // println!("Liquid Volume Fraction (-) = {:.3}", p1.RL);
    // println!("Two-Phase Velocity (m/sec) = {:.4}", p1.UTP);
    // println!("1.0 Velocity Head (kgf/cm^2) = {:.4}", p1.Head);
    // println!("Frictional Pressure Loss (kgf/cm^2/100m) = {:.4}", p1.Pfric);
    // println!("Elevation Head Loss (kgf/cm^2/100m) = {:.4}", p1.Pgrav);
    // println!("Erosion Factor (-) = {:.3}", p1.Ef);
    // println!("if Φ ≤ 1 : No Erosion, Φ > 1, Erosion occurred");
    // EndRegion

    //Region Test data for Stratified Wavy Flow (Stratified Model)
    // Liquid data
    // let wl: f64 = 64870.62744; // [kg/hr]
    // let lo_l: f64 = 790.9917; // [kg/m^3]
    // let mu_l: f64 = 0.241; // [cP]
    // let surface_tension: f64 = 14.78; // [dyne/cm]
    //                                   // Vapor data
    // let wg: f64 = 21623.54248; // [kg/hr]
    // let lo_g: f64 = 4.58128; // [kg/m^3]
    // let mu_g: f64 = 0.0091; // [cP]
    //                         // Misc. data
    // let id: f64 = 23.25; // [in]
    // let slope: f64 = 0.0; // [degree]
    // let rough: f64 = 0.04572; // [mm]
    // let sf: f64 = 1.0; // [-]
    //
    // let mut p2 = Horizontal::new(
    //     wl,
    //     wg,
    //     lo_l,
    //     lo_g,
    //     mu_l,
    //     mu_g,
    //     surface_tension,
    //     rough,
    //     sf,
    //     id,
    //     slope,
    // );
    // p2.unit_transfer();
    // p2.flow_regime();
    // println!("p2 flow regime << {} >>", p2.flow_regime);
    // p2.model_cal();
    // println!("Tow-Phase Density (kg/m^3) = {:.4}", p2.LoTP);
    // println!("Liquid Depth-BOP (m) = {:.4}", p2.depth);
    // println!("Liquid Velocity (m/sec) = {:.4}", p2.velL);
    // println!("Vapor Velocity (m/sec) = {:.4}", p2.velG);
    // println!("Liquid Volume Fraction (-) = {:.4}", p2.RL);
    // println!("1.0 Velocity Head (kgf/cm^2) = {:.4}", p2.Head);
    // println!("Frictional Pressure Loss (kgf/cm^2/100m) = {:.4}", p2.Pfric);
    // println!("Erosion Factor (-) = {:.3}", p2.Ef);
    // println!("if Φ ≤ 1 : No Erosion, Φ > 1, Erosion occurred");
    // EndRegion

    //Region Test data for Intermittent-Slug Flow (Slug Model)
    // Liquid data
    // let wl: f64 = 116604.0; // [kg/hr]
    // let lo_l: f64 = 803.0; // [kg/m^3]
    // let mu_l: f64 = 0.45; // [cP]
    // let surface_tension: f64 = 9.3; // [dyne/cm]
    //                                 // Vapor data
    // let wg: f64 = 1896.0; // [kg/hr]
    // let lo_g: f64 = 8.49; // [kg/m^3]
    // let mu_g: f64 = 0.02; // [cP]
    //                       // Misc. data
    // let id: f64 = 7.981; // [in]
    // let slope: f64 = 0.0; // [degree]
    // let rough: f64 = 0.046; // [mm]
    // let sf: f64 = 1.0; // [-]
    //
    // let mut p3 = Horizontal::new(
    //     wl,
    //     wg,
    //     lo_l,
    //     lo_g,
    //     mu_l,
    //     mu_g,
    //     surface_tension,
    //     rough,
    //     sf,
    //     id,
    //     slope,
    // );
    // p3.unit_transfer();
    // p3.flow_regime();
    // println!("p3 flow regime << {} >>", p3.flow_regime);
    // p3.model_cal();
    // println!("Tow-Phase Slug Unit Density (kg/m^3) = {:.4}", p3.LoSU);
    // println!("Liquid Slug Unit Density [m/s] = {:.4}", p3.LoLS);
    // println!("Liquid Volume Fraction (-) = {:.4}", p3.RL);
    // println!("Liquid Slug Velocity (m/sec) = {:.4}", p3.Us);
    // println!("Liquid Slug Length (m) = {:.4}", p3.Ls);
    // println!("Liquid Slug Length (m) = {:.4}", p3.Lu);
    // println!("1.0 Velocity Head (kgf/cm^2) = {:.4}", p3.Head);
    // println!("Frictional Pressure Loss (kgf/cm^2/100m) = {:.4}", p3.Pfric);
    // println!("Erosion Factor (-) = {:.3}", p3.Ef);
    // println!("if Φ ≤ 1 : No Erosion, Φ > 1, Erosion occurred");
    // EndRegion
}

pub fn vertical_down_validate() {
    //Region Test data for Annular Flow
    // Liquid data
    let wl: f64 = 21937.88; // [kg/hr]
    let lo_l: f64 = 962.0689; // [kg/m^3]
    let mu_l: f64 = 0.511; // [cP]
    let surface_tension: f64 = 30.3; // [dyne/cm]
                                     // Vapor data
    let wg: f64 = 376.93329; // [kg/hr]
    let lo_g: f64 = 0.9931447; // [kg/m^3]
    let mu_g: f64 = 0.01; // [cP]
                          // Misc. data
    let id: f64 = 15.25; // [in]
    let slope: f64 = 0.0; // [degree]
    let rough: f64 = 0.04572; // [mm]
    let sf: f64 = 1.0; // [-]

    let mut p1 = VerticalDown::new(
        wl,
        wg,
        lo_l,
        lo_g,
        mu_l,
        mu_g,
        surface_tension,
        rough,
        sf,
        id,
        slope,
    );
    p1.unit_transfer();
    p1.flow_regime();
    println!("flow regime << {} >>", p1.flow_regime);
    // p1.model_cal();
    // println!("Tow-Phase Density (kg/m^3) = {:.4}", p1.Loip);
    // println!("Liquid Volume Fraction (-) = {:.3}", p1.RL);
    // println!("Two-Phase Velocity (m/sec) = {:.4}", p1.UTP);
    // println!("1.0 Velocity Head (kgf/cm^2) = {:.4}", p1.Head);
    // println!("Frictional Pressure Loss (kgf/cm^2/100m) = {:.4}", p1.Pfric);
    // println!("Elevation Head Loss (kgf/cm^2/100m) = {:.4}", p1.Pgrav);
    // println!("Erosion Factor (-) = {:.3}", p1.Ef);
    // println!("if Φ ≤ 1 : No Erosion, Φ > 1, Erosion occurred");
    // EndRegion

    //Region Test data for Dispersed Bubble Flow
    // Liquid data
    // let wl: f64 = 4536.0; // [kg/hr]
    // let lo_l: f64 = 801.0; // [kg/m^3]
    // let mu_l: f64 = 0.6; // [cP]
    // let surface_tension: f64 = 10.0; // [dyne/cm]
    //                                  // Vapor data
    // let wg: f64 = 0.4536; // [kg/hr]
    // let lo_g: f64 = 8.0; // [kg/m^3]
    // let mu_g: f64 = 0.01; // [cP]
    //                       // Misc. data
    // let id: f64 = 1.049; // [in]
    // let slope: f64 = 0.0; // [degree]
    // let rough: f64 = 0.046; // [mm]
    // let sf: f64 = 1.0; // [-]
    //
    // let mut p2 = VerticalDown::new(
    //     wl,
    //     wg,
    //     lo_l,
    //     lo_g,
    //     mu_l,
    //     mu_g,
    //     surface_tension,
    //     rough,
    //     sf,
    //     id,
    //     slope,
    // );
    // p2.unit_transfer();
    // p2.flow_regime();
    // println!("flow regime << {} >>", p2.flow_regime);
    // p2.model_cal();
    // println!("Tow-Phase Density (kg/m^3) = {:.4}", p2.Loip);
    // println!("Liquid Volume Fraction (-) = {:.3}", p2.RL);
    // println!("Two-Phase Velocity (m/sec) = {:.4}", p2.UTP);
    // println!("1.0 Velocity Head (kgf/cm^2) = {:.4}", p2.Head);
    // println!("Frictional Pressure Loss (kgf/cm^2/100m) = {:.4}", p2.Pfric);
    // println!("Elevation Head Loss (kgf/cm^2/100m) = {:.4}", p2.Pgrav);
    // println!("Erosion Factor (-) = {:.3}", p2.Ef);
    // println!("if Φ ≤ 1 : No Erosion, Φ > 1, Erosion occurred");
    // EndRegion

    //Region Test data for Slug Flow
    // Liquid data
    // let wl: f64 = 12000.0; // [kg/hr]
    // let lo_l: f64 = 801.0; // [kg/m^3]
    // let mu_l: f64 = 0.6; // [cP]
    // let surface_tension: f64 = 10.0; // [dyne/cm]
    //                                  // Vapor data
    // let wg: f64 = 450.0; // [kg/hr]
    // let lo_g: f64 = 8.0; // [kg/m^3]
    // let mu_g: f64 = 0.01; // [cP]
    //                       // Misc. data
    // let id: f64 = 1.049; // [in]
    // let slope: f64 = 0.0; // [degree]
    // let rough: f64 = 0.046; // [mm]
    // let sf: f64 = 1.0; // [-]
    //
    // let mut p3 = VerticalDown::new(
    //     wl,
    //     wg,
    //     lo_l,
    //     lo_g,
    //     mu_l,
    //     mu_g,
    //     surface_tension,
    //     rough,
    //     sf,
    //     id,
    //     slope,
    // );
    // p3.unit_transfer();
    // p3.flow_regime();
    // println!("flow regime << {} >>", p3.flow_regime);
    // p1.model_cal();
    // println!("Tow-Phase Density (kg/m^3) = {:.4}", p3.Loip);
    // println!("Liquid Volume Fraction (-) = {:.3}", p3.RL);
    // println!("Two-Phase Velocity (m/sec) = {:.4}", p3.UTP);
    // println!("1.0 Velocity Head (kgf/cm^2) = {:.4}", p3.Head);
    // println!("Frictional Pressure Loss (kgf/cm^2/100m) = {:.4}", p3.Pfric);
    // println!("Elevation Head Loss (kgf/cm^2/100m) = {:.4}", p3.Pgrav);
    // println!("Erosion Factor (-) = {:.3}", p3.Ef);
    // println!("if Φ ≤ 1 : No Erosion, Φ > 1, Erosion occurred");
    // EndRegion
}

fn main() {
    // vertical_up_validate();
    // horizontal_validate();
    vertical_down_validate();
}
