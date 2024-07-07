#[derive(Debug)]
pub enum Regime {
    // Vertical Up Flow Regime
    VerticalUpAnnularFlow(String),
    VerticalUpBubbleFlow(String),
    VerticalUpSlugAndChurnFlow(String),
    VerticalUpFinelyDispersedBubbleFlow(String),
    // Horizontal Flow Regime
    HorizontalStratifiedSmoothFlow(String),
    HorizontalStratifiedWavyFlow(String),
    HorizontalAnnularDispersedFlow(String),
    HorizontalElongatedBubbleFlow(String),
    HorizontalIntermittentSlugFlow(String),
    HorizontalDispersedBubbleFlow(String),
    // Vertical Down Flow Regime
    VerticalDownAnnularFlow(String),
    VerticalDownSlugFlow(String),
    VerticalDownDispersedBubbleFlow(String),
    // Others
    NONE,
}
pub trait TwoPhaseLine {
    fn unit_transfer(&mut self);
    fn flow_regime(&mut self);

    fn model_cal(&mut self);
}
