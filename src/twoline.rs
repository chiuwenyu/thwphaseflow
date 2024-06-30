#[derive(Debug)]
pub enum Regime {
    VerticalUpAnnularFlow(String),
    VerticalUpBubbleFlow(String),
    VerticalUpSlugAndChurnFlow(String),
    VerticalUpFinelyDispersedBubbleFlow(String),
    NONE,
}
pub trait TwoPhaseLine {
    fn unit_transfer(&mut self);
    fn flow_regime(&mut self);
    #[allow(dead_code)]
    fn model_cal();
}