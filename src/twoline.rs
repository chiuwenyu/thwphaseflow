#[derive(Debug)]
pub enum Regime {
    VerticalUpAnnularFlow(String),
    VerticalUpBubbleFlow(String),
    VerticalUpSlugAndChurnFlow(String),
    VerticalUpFinelyDispersedBubbleFlow(String),
    HorizontalStratifiedSmoothFlow(String),
    HorizontalStratifiedWavyFlow(String),
    HorizontalAnnularDispersedFlow(String),
    HorizontalElongatedBubbleFlow(String),
    HorizontalIntermittentSlugFlow(String),
    HorizontalDispersedBubbleFlow(String),
    NONE,
}
pub trait TwoPhaseLine {
    fn unit_transfer(&mut self);
    fn flow_regime(&mut self);

    fn model_cal(&mut self);
}