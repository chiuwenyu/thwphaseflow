pub enum Regime {
    VerticalUpAnnularFlow(String),
    VerticalUpBubbleFlow(String),
    VerticalUpSlugAndChurnFlow(String),
    VerticalUpFinelyDispersedBubbleFlow(String),
    NONE,
}
pub trait TwoPhaseLine {
    fn unit_transfer(&self);
    fn flow_regime(&self) -> Regime;
    fn model_cal();
}