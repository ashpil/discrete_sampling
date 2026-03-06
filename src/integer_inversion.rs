use crate::distribution::{
    Discrete1D,
    Discrete1DPdf,
    Continuous1D,
};
use num_traits::{
    Num,
    real::Real,
    AsPrimitive,
};
#[cfg(not(feature = "std"))]
use alloc::boxed::Box;

use crate::Inversion1D;

pub type IntegerInversion2D<R> = crate::Adapter2D<IntegerInversion1D<R>>;

#[cfg_attr(feature = "rkyv", derive(rkyv::Archive, rkyv::Deserialize, rkyv::Serialize))]
pub struct IntegerInversion1D<W> {
    base: Inversion1D<W>,
    inverse_cdf: Box<[usize]>,
}

impl<W: Num + AsPrimitive<R> + PartialOrd, R: Real + AsPrimitive<W> + AsPrimitive<usize> + 'static> Discrete1D<R> for IntegerInversion1D<W>
    where usize: AsPrimitive<W>,
          usize: AsPrimitive<R>,
{
    type Weight = W;

    fn build(weights: &[W]) -> Self {
        let base = Inversion1D::build(weights);

        let mut lower_bound = 0;
        let size = base.cdf.len(); // technically this is a free parameter, but at this size the algorithm becomes expected O(1) for uniform inputs
        let inverse_cdf = (0..size).map(|i| {
            while base.cdf[lower_bound + 1] * size.as_() < <usize as AsPrimitive<W>>::as_(i) * *base.cdf.last().unwrap() {
                lower_bound += 1;
            }
            lower_bound
        }).collect::<Box<[usize]>>();

        Self {
            base,
            inverse_cdf,
        }
    }

    fn sample(&self, u: R) -> usize {
        let point: W = (u * self.integral().as_()).as_();
        let start_offset: usize = (u * self.inverse_cdf.len().as_()).as_();
        let mut index = self.inverse_cdf[start_offset];
        while self.base.cdf[index + 1] <= point {
            index += 1;
        }
        index
    }

    fn integral(&self) -> W {
        self.base.integral()
    }

    fn size(&self) -> usize {
        self.base.size()
    }
}

impl<W: Num + AsPrimitive<R> + PartialOrd, R: Real + AsPrimitive<W> + AsPrimitive<usize> + 'static> Discrete1DPdf<R> for IntegerInversion1D<W>
    where usize: AsPrimitive<W>,
          usize: AsPrimitive<R>,
{
    fn pdf(&self, u: usize) -> W {
        self.base.pdf(u)
    }
}

impl<W: Num + AsPrimitive<R> + PartialOrd, R: Real + AsPrimitive<W> + AsPrimitive<usize> + 'static> Continuous1D<R> for IntegerInversion1D<W>
    where usize: AsPrimitive<W>,
          usize: AsPrimitive<R>,
{
    fn sample_continuous(&self, u: R) -> R {
        let offset = self.sample(u);
        let du = (u * self.integral().as_() - self.base.cdf[offset].as_()) / (self.base.cdf[offset + 1].as_() - self.base.cdf[offset].as_());
        (<usize as AsPrimitive<R>>::as_(offset) + du) / self.size().as_()
    }

    fn invert_continuous(&self, u: R) -> R {
        let scaled: R = <usize as AsPrimitive<R>>::as_(self.size()) * u;
        let idx: usize = scaled.as_();
        let delta = scaled - idx.as_();
        crate::utils::lerp(delta, self.base.cdf[idx].as_(), self.base.cdf[idx + 1].as_()) / self.integral().as_()
    }
}

#[cfg(test)]
mod tests {
    use crate::distribution::distribution_1d_tests;
    use crate::distribution::continuous_distribution_1d_tests;

    distribution_1d_tests!(crate::integer_inversion::IntegerInversion1D);
    continuous_distribution_1d_tests!(crate::integer_inversion::IntegerInversion1D);
}

