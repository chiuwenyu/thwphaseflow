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

    fn model_cal(&self);
}