#[macro_export]
macro_rules! sin_deg {
    ($x: expr) => {
        $x.to_radians().sin()
    };
}

#[macro_export]
macro_rules! cos_deg {
    ($x: expr) => {
        $x.to_radians().cos()
    };
}

#[macro_export]
macro_rules! atan2_deg {
    ($x : expr, $y : expr) => {
        $x.atan2($y).to_degrees()
    };
}

#[macro_export]
macro_rules! in_range {
    ($lb: expr, $ub: expr, $val: expr) => {
        $val >= $lb && $val <= $ub
    };
}
