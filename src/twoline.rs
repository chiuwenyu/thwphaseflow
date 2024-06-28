pub enum Regime {
    VerticalUpAnnularFlow,
    VerticalUpChurnFlow,
    VerticalUpBubbleFlow,
    VerticalUpSlugAndChurnFlow,
    VerticalUpFinelyDispersedBubbleFlow,
    NONE,
}
pub trait TwoPhaseLine {
    fn unit_transfer(&self);
    fn flow_regime(&self) -> Regime;
    fn model_cal();
}