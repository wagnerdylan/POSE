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
    ($y : expr, $x : expr) => {
        $y.atan2($x).to_degrees()
    };
}

#[macro_export]
macro_rules! in_range {
    ($lb: expr, $ub: expr, $val: expr) => {
        $val >= $lb && $val <= $ub
    };
}
