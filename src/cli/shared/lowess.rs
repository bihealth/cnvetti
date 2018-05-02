use std::cmp;

pub fn pow2(x: f64) -> f64 {
    x * x
}

pub fn pow3(x: f64) -> f64 {
    x * x * x
}

#[derive(Debug)]
pub struct LowessResults {
    // Fitted value.
    pub fit: Vec<f64>,
    // Residuals.
    pub residuals: Vec<f64>,
    // Residual weights.
    pub residual_weights: Vec<f64>,
}

// Internal helper structure for helping in the reversible sorting.
struct LoessElem {
    x: f64,
    y: f64,
    no: usize,
}

/// Rust implementation of `LOWESS()` based on RATFOR code translated from original FORTRAN
/// code.
///
/// # Args
///
/// - `x` - The x coordinates in the scatter plot.
/// - `y` - The y coordinates in the scatter plot.
/// - `f` - Smoother span, R default is 2/3.
/// - `nsteps` - Number of iteration steps, R default is 3.
/// - `delta` - The `delta` value, R default is to use 1/100th of the range of `x`.
///
/// # Result
///
/// Struct `LowessResults` reflecting the fit, residuals, and residual weights.
pub fn lowess(x: &[f64], y: &[f64], f: f64, nsteps: usize, delta: f64) -> LowessResults {
    assert_eq!(x.len(), y.len());

    let mut elems: Vec<LoessElem> = x.iter()
        .zip(y.iter())
        .enumerate()
        .map(|(no, (x, y))| LoessElem {
            x: *x,
            y: *y,
            no: no,
        })
        .collect();
    elems.sort_by(|a, b| a.x.partial_cmp(&b.x).unwrap());
    let x: Vec<f64> = elems.iter().map(|elem| elem.x).collect();
    let y: Vec<f64> = elems.iter().map(|elem| elem.y).collect();

    let n = x.len();
    let mut result = LowessResults {
        fit: vec![y[0]; n],
        residual_weights: vec![0.0; n],
        residuals: vec![0.0; n],
    };

    if n < 2 {
        return result;
    }

    // At least two, at most n points.
    let ns = cmp::max(2, cmp::min(((f * n as f64) + 1e-7) as usize, n));
    // Robustness iteration.
    for iter in 1..(nsteps + 1) {
        println!("iter = {}", iter);
        let mut nleft = 0;
        let mut nright = ns - 1;
        let mut last = 0; // index of prev estimated point
        let mut i = 1; // index of current point
        loop {
            if nright < n - 1 {
                // Move nleft, nright to right if radius decreases.
                let d1 = x[i] - x[nleft];
                let d2 = x[nright + 1] - x[i];
                if d2 > d1 {
                    // Radius will not decrease by move right.
                    nleft += 1;
                    nright += 1;
                    continue;
                }
            }

            // Compute fitted value at x[i].
            let l = lowest(
                &x,
                &y,
                x[i],
                nleft,
                nright,
                &mut result.residuals,
                if iter > 1 {
                    Some(&mut result.residual_weights)
                } else {
                    None
                },
            );
            result.fit[i] = match l {
                Ok(val) => val,
                Err(_) => y[i],
            };

            // All weights zero - copy over value (all rw==0)
            if last < i - 1 {
                // Skipped points -- interpolate.
                let denom = x[i] - x[last];
                for j in (last + 1)..i {
                    let alpha = (x[j] - x[last]) / denom;
                    result.fit[j] = alpha * result.fit[i] + (1.0 - alpha) * result.fit[last];
                }
            }
            // Last point actually estimated.
            last = i;
            // x coord of close points.
            let cut = x[i] + delta;
            // Find close points.
            for i in (last + 1)..n {
                if x[i] > cut {
                    break; // i one beyond the last point within cut
                } else if x[i] == x[last] {
                    // Exact match in x.
                    result.fit[i] = result.fit[last];
                    last = i;
                }
            }

            // Back 1 point so interpolation within delta, but always go forward.
            i = cmp::max(last + 1, i - 1);

            if last >= n - 1 {
                break;
            }
        }

        for i in 0..n {
            // Residuals.
            result.residuals[i] = y[i] - result.fit[i];
        }

        // Overall scal estimate
        let mut sc = 0.0;
        for i in 0..n {
            sc += result.residuals[i].abs();
        }
        sc = sc / (n as f64);

        // Compute robustness weights excapt last time.
        if iter > nsteps {
            break;
        }

        for i in 0..n {
            result.residual_weights[i] = result.residuals[i].abs();
        }

        result
            .residual_weights
            .sort_by(|a, b| a.partial_cmp(b).unwrap());

        let m1 = 1 + n / 2;
        let m2 = n - m1 + 1;
        // Compute 6 median abs residual.
        let cmad = 3.0 * (result.residual_weights[m1] + result.residual_weights[m2]);
        if cmad < 1e-7 * sc {
            break; // effectively zero, trick from R implementation
        }
        let c9 = 0.999 * cmad;
        let c1 = 0.001 * cmad;
        for i in 0..n {
            let r = result.residuals[i].abs();
            if r <= c1 {
                result.residual_weights[i] = 1.0; // near 0, avoid underflow
            } else if r > c9 {
                result.residual_weights[i] = 0.0; // near 1, avoid underflow
            } else {
                result.residual_weights[i] = pow2(1.0 - pow2(r / cmad));
            }
        }
    }

    // Permute back the values into the original order.
    let mut tmp = LowessResults {
        fit: vec![y[0]; n],
        residual_weights: vec![0.0; n],
        residuals: vec![0.0; n],
    };
    for (i, no) in elems.iter().map(|elem| elem.no).enumerate() {
        tmp.fit[no] = result.fit[i];
        tmp.residual_weights[no] = result.residual_weights[i];
        tmp.residuals[no] = result.residuals[i];
    }

    tmp
}

pub fn lowest(
    x: &[f64],
    y: &[f64],
    xs: f64,
    nleft: usize,
    nright: usize,
    w: &mut [f64],
    rw: Option<&[f64]>,
) -> Result<f64, ()> {
    let n = x.len();
    let range = x.last().unwrap() - x.first().unwrap();
    let h = (xs - x[nleft]).max(x[nright] - xs);
    let h9 = 0.999 * h;
    let h1 = 0.001 * h;

    // Compute weights (pick up all ties on right)
    let mut a = 0.0_f64;
    let mut j = nleft;
    while j < n {
        w[j] = 0.0;
        let r = (x[j] - xs).abs();
        if r <= h9 {
            // small enough for non-zero weight
            if r > h1 {
                w[j] = pow3(1.0 - pow3(r / h));
            } else {
                w[j] = 1.0;
            }
            if let Some(ref rw) = rw {
                w[j] = rw[j] * w[j];
            }
            a += w[j];
        } else if x[j] > xs {
            break; // get out at first zero weight on right
        }
        j += 1;
    }

    let nrt = j - 1; // rightmost ponit (may be greater than nright because of ties)
    if a <= 0.0 {
        return Err(());
    }

    // Make sum of w[j] == 1.
    for j in nleft..(nrt + 1) {
        w[j] = w[j] / a;
    }
    if h > 0.0 {
        // use linear fit
        // Find weighted center of x values
        let mut a = 0.0;
        for j in nleft..(nrt + 1) {
            a += w[j] * x[j];
        }

        let mut b = xs - a;
        let mut c = 0.0_f64;
        for j in nleft..(nrt + 1) {
            c += w[j] * pow2(x[j] - a);
        }

        if c.sqrt() > 0.001 * range {
            // Points are spread out enough to compute slope.
            b = b / c;
            for j in nleft..(nrt + 1) {
                w[j] = w[j] * (1.0 + b * (x[j] - a));
            }
        }
    }

    let mut ys = 0.0;
    for j in nleft..(nrt + 1) {
        ys += w[j] * y[j];
    }

    Ok(ys)
}

/*
    The lowess code above is Translated from RATFOR lowess code of W. S.
    Cleveland as obtained from NETLIB.

    It is based on two functions written in ratfor (see below), namely lowest
    and lowess. The code has since been refactored and commented further.
  */

/* ratfor code for lowest:
  *
  *  subroutine lowest(x,y,n,xs,ys,nleft,nright,w,userw,rw,ok)
  *  real x(n),y(n),w(n),rw(n)
  *  logical userw,ok
  *  range = x(n)-x(1)
  *  h = amax1(xs-x(nleft),x(nright)-xs)
  *  h9 = .999*h
  *  h1 = .001*h
  *  a = 0.0        # sum of weights
  *  for(j=nleft; j<=n; j=j+1){     # compute weights (pick up all ties on right)
  *         w(j)=0.
  *         r = abs(x(j)-xs)
  *         if (r<=h9) {    # small enough for non-zero weight
  *                 if (r>h1) w(j) = (1.0-(r/h)**3)**3
  *                 else      w(j) = 1.
  *                 if (userw) w(j) = rw(j)*w(j)
  *                 a = a+w(j)
  *                 }
  *         else if(x(j)>xs)break   # get out at first zero wt on right
  *         }
  *  nrt=j-1        # rightmost pt (may be greater than nright because of ties)
  *  if (a<=0.0) ok = FALSE
  *  else { # weighted least squares
  *         ok = TRUE
  *         do j = nleft,nrt
  *                 w(j) = w(j)/a   # make sum of w(j) == 1
  *         if (h>0.) {     # use linear fit
  *                 a = 0.0
  *                 do j = nleft,nrt
  *                         a = a+w(j)*x(j) # weighted center of x values
  *                 b = xs-a
  *                 c = 0.0
  *                 do j = nleft,nrt
  *                         c = c+w(j)*(x(j)-a)**2
  *                 if(sqrt(c)>.001*range) {
  *  # points are spread out enough to compute slope
  *                         b = b/c
  *                         do j = nleft,nrt
  *                                 w(j) = w(j)*(1.0+b*(x(j)-a))
  *                         }
  *                 }
  *         ys = 0.0
  *         do j = nleft,nrt
  *                 ys = ys+w(j)*y(j)
  *         }
  *  return
  *  end

  */

/* ratfor code for lowess:
 *
 *  subroutine lowess(x,y,n,f,nsteps,delta,ys,rw,res)
 *  real x(n),y(n),ys(n),rw(n),res(n)
 *  logical ok
 *  if (n<2){ ys(1) = y(1); return }
 *  ns = max0(min0(ifix(f*float(n)),n),2)  # at least two, at most n points
 *  for(iter=1; iter<=nsteps+1; iter=iter+1){      # robustness iterations
 *         nleft = 1; nright = ns
 *         last = 0        # index of prev estimated point
 *         i = 1   # index of current point
 *         repeat{
 *                 while(nright<n){
 *  # move nleft, nright to right if radius decreases
 *                         d1 = x(i)-x(nleft)
 *                         d2 = x(nright+1)-x(i)
 *  # if d1<=d2 with x(nright+1)==x(nright), lowest fixes
 *                         if (d1<=d2) break
 *  # radius will not decrease by move right
 *                         nleft = nleft+1
 *                         nright = nright+1
 *                         }
 *                 call lowest(x,y,n,x(i),ys(i),nleft,nright,res,iter>1,rw,ok)
 *  # fitted value at x(i)
 *                 if (!ok) ys(i) = y(i)
 *  # all weights zero - copy over value (all rw==0)
 *                 if (last<i-1) { # skipped points -- interpolate
 *                         denom = x(i)-x(last)    # non-zero - proof?
 *                         for(j=last+1; j<i; j=j+1){
 *                                 alpha = (x(j)-x(last))/denom
 *                                 ys(j) = alpha*ys(i)+(1.0-alpha)*ys(last)
 *                                 }
 *                         }
 *                 last = i        # last point actually estimated
 *                 cut = x(last)+delta     # x coord of close points
 *                 for(i=last+1; i<=n; i=i+1){     # find close points
 *                         if (x(i)>cut) break     # i one beyond last pt within cut
 *                         if(x(i)==x(last)){      # exact match in x
 *                                 ys(i) = ys(last)
 *                                 last = i
 *                                 }
 *                         }
 *                 i=max0(last+1,i-1)
 *  # back 1 point so interpolation within delta, but always go forward
 *                 } until(last>=n)
 *         do i = 1,n      # residuals
 *                 res(i) = y(i)-ys(i)
 *         if (iter>nsteps) break  # compute robustness weights except last time
 *         do i = 1,n
 *                 rw(i) = abs(res(i))
 *         call sort(rw,n)
 *         m1 = 1+n/2; m2 = n-m1+1
 *         cmad = 3.0*(rw(m1)+rw(m2))      # 6 median abs resid
 *         c9 = .999*cmad; c1 = .001*cmad
 *         do i = 1,n {
 *                 r = abs(res(i))
 *                 if(r<=c1) rw(i)=1.      # near 0, avoid underflow
 *                 else if(r>c9) rw(i)=0.  # near 1, avoid underflow
 *                 else rw(i) = (1.0-(r/cmad)**2)**2
 *                 }
 *         }
 *  return
 *  end
 */

#[cfg(test)]
mod tests {
    use super::lowess;

    #[test]
    fn test_lowest_cars() {
        let speed = vec![
            4, 4, 7, 7, 8, 9, 10, 10, 10, 11, 11, 12, 12, 12, 12, 13, 13, 13, 13, 14, 14, 14, 14,
            15, 15, 15, 16, 16, 17, 17, 17, 18, 18, 18, 18, 19, 19, 19, 20, 20, 20, 20, 20, 22, 23,
            24, 24, 24, 24, 25,
        ];
        let fspeed: Vec<f64> = speed.iter().map(|x| *x as f64).collect();
        let dist = vec![
            2, 10, 4, 22, 16, 10, 18, 26, 34, 17, 28, 14, 20, 24, 28, 26, 34, 34, 46, 26, 36, 60,
            80, 20, 26, 54, 32, 40, 32, 40, 50, 42, 56, 76, 84, 36, 46, 68, 32, 48, 52, 56, 64, 66,
            54, 70, 92, 93, 120, 85,
        ];
        let fdist: Vec<f64> = dist.iter().map(|x| *x as f64).collect();

        let res = lowess(
            &fspeed,
            &fdist,
            0.9,
            3,
            0.01 * ((dist.iter().max().unwrap() - dist.iter().min().unwrap()) as f64),
        );

        println!("{:?}", res);
    }
}
