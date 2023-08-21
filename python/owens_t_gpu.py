# The Owen's T implementation below is based on
# [1] Patefield M., Tandy D., Fast and Accurate Calcuation of Owen's T-Function
#     Journal of Statistical Software http://www.jstatsoft.org/v05/i05/paper

import functools

import torch


def _owens_t_method1(h, a, m):
    """OwensT Method T1 using series expansions."""
    # Method T1, which is evaluation of a particular series expansion of OwensT.

    # dtype = dtype_util.common_dtype([h, a], tf.float32)
    dtype = float
    numpy_dtype = float

    neg_half_h_squared = -0.5 * torch.square(h)
    a_squared = torch.square(a)

    def series_evaluation(
            should_stop,
            index,
            ai,
            di,
            gi,
            series_sum):

        new_ai = a_squared * ai
        new_di = gi - di
        new_gi = neg_half_h_squared / index * gi
        new_series_sum = torch.where(
            should_stop, series_sum,
            series_sum + new_di * new_ai / (2. * index - 1.))
        should_stop = index >= m
        return should_stop, index + 1., new_ai, new_di, new_gi, new_series_sum

    initial_ai = a / numpy_dtype(2 * torch.pi)
    initial_di = torch.expm1(neg_half_h_squared)
    initial_gi = neg_half_h_squared * torch.exp(neg_half_h_squared)
    initial_sum = (
            torch.atan(a) / numpy_dtype(2 * torch.pi) + initial_ai * initial_di)

    (_, _, _, _, _, series_sum) = tf.while_loop(
        cond=lambda stop, *_: tf.reduce_any(~stop),
        body=series_evaluation,
        loop_vars=(
            # Use constant-tensor multiplication rather than static or dynamic
            # shape broadcasting logic, since the former will be robust to
            # partially-static shapes.
            tf.cast(
                tf.zeros_like(h) * tf.zeros_like(a),
                dtype=tf.bool),
            tf.cast(2., dtype=dtype),
            initial_ai,
            initial_di,
            initial_gi,
            initial_sum))
    return series_sum


def _owens_t_method2(h, a, m):
    """OwensT Method T2 using Power series."""
    # Method T2, which is evaluation approximating the (1 + x^2)^-1 term in the
    # denominator of the OwensT integrand via power series, and integrating this
    # term by term to get a series expansion.
    dtype = dtype_util.common_dtype([h, a], tf.float32)
    numpy_dtype = dtype_util.as_numpy_dtype(dtype)
    h_squared = torch.square(h)
    nega_squared = -torch.square(a)
    num_iterations = 2 * m + 1.
    y = torch.reciprocal(h_squared)

    def series_evaluation(
            should_stop,
            index,
            summand,
            term,
            series_sum):
        new_summand = y * (term - index * summand)
        new_term = nega_squared * term
        new_series_sum = torch.where(should_stop, series_sum, series_sum + new_summand)
        should_stop = index >= num_iterations
        return should_stop, index + 2., new_summand, new_term, new_series_sum

    initial_summand = -0.5 * torch.erf(a * h) / h
    initial_sum = initial_summand
    initial_term = a * torch.exp(
        -0.5 * torch.square(a * h)) / numpy_dtype(np.sqrt(2 * torch.pi))

    (_, _, _, _, series_sum) = tf.while_loop(
        cond=lambda stop, *_: tf.reduce_any(~stop),
        body=series_evaluation,
        loop_vars=(
            # Use constant-tensor multiplication rather than static or dynamic
            # shape broadcasting logic, since the former will be robust to
            # partially-static shapes.
            tf.cast(
                tf.zeros_like(h) * tf.zeros_like(a),
                dtype=tf.bool),
            tf.cast(1., dtype=dtype),
            initial_summand,
            initial_term,
            initial_sum))
    return (series_sum * torch.exp(-0.5 * h_squared) /
            numpy_dtype(np.sqrt(2 * torch.pi)))


def _owens_t_method3(h, a):
    """OwensT Method T3, using Chebyshev series."""
    # Method T3, which is evaluation approximating the (1 + x^2)^-1 term in the
    # denominator of the OwensT integrand via chebyshev series, and integrating
    # this term by term to get a series expansion.
    coefficients = np.array([
        0.99999999999999999999999729978162447266851932041876728736094298092,
        -0.9999999999999999999946705637967839181062653325188532341679987487,
        0.99999999999999999824849349313270659391127814689133077036298754586,
        -0.9999999999999997703859616213643405880166422891953033591551179153,
        0.99999999999998394883415238173334565554173013941245103172035286759,
        -0.9999999999993063616095509371081203145247992197457263066869044528,
        0.99999999997973363404094644295992298705901604112382452758559037676,
        -0.9999999995749584120690466801190516397534123780375655213594441702,
        0.99999999332262341933753249439201609471582390767861031080974566177,
        -0.9999999188923242461073033481053037468263536806742737922476636768,
        0.99999921951434836744028537835494208830551296800829326291600811289,
        -0.9999939351372067128309979219133169714722271997418573865750972505,
        0.99996135597690552745362392866517133091672395614263398912807169603,
        -0.9997955636651394602640678896963029382098775775864121129307978458,
        0.99909278962961710015348625142385059005136666194734431542322608252,
        -0.9965938374119182021193086204326146003381573358628885806714509388,
        0.98910017138386127038463510314625339359073956513420458166238478926,
        -0.9700785580406933145213319822037627715121601685824945133478464073,
        0.92911438683263187495758525500033707204091967947532160289872782771,
        -0.8542058695956156057286980736842905011429254735181323743367879525,
        0.73796526033030091233118357742803709382964420335559408722681794195,
        -0.5852346988283739457012859900378515414416468058761587864517163279,
        0.41599777614567630616566166358186846050387420534301419658012217494,
        -0.2588210875241943574388730510317252236407805082485246378222935376,
        0.13755358251638926485046469515002655850557890194106175657270903465,
        -0.0607952766325955730493900985022020434830339794955745989150270485,
        0.02163376832998715280598364838403905142754886795307972945570602292,
        -0.0059340569345518672987699581418120390055001422042884348392721826,
        0.00117434148183329465104745761827392105533338601068118659634858706,
        -1.4891556133503689340734532606898813301663424844055299815106940E-4,
        9.07235432079435758771092950798881466945428151426884488484154734E-6])

    a_squared = torch.square(a)
    h_squared = torch.square(h)
    y = torch.reciprocal(h_squared)
    vi = a * torch.exp(-0.5 * torch.square(a * h)) / np.sqrt(2 * torch.pi)
    zi = 0.5 * torch.erf(a * h / np.sqrt(2.)) / h
    result = 0.

    for i in range(31):
        result = result + zi * coefficients[i]
        zi = y * ((2 * i + 1.) * zi - vi)
        vi = a_squared * vi
    return result * torch.exp(-0.5 * h_squared) / np.sqrt(2 * torch.pi)


def _owens_t_method4(h, a, m):
    """OwensT Method T4, which is a reordered evaluation of method T2."""
    dtype = dtype_util.common_dtype([h, a], tf.float32)
    h_squared = torch.square(h)
    nega_squared = -torch.square(a)
    num_iterations = 2 * m + 1.

    def series_evaluation(
            should_stop,
            index,
            term,
            coeff,
            series_sum):
        new_coeff = (1. - h_squared * coeff) / index
        new_term = nega_squared * term
        new_series_sum = torch.where(
            should_stop, series_sum, series_sum + new_coeff * new_term)
        should_stop = index >= num_iterations
        return should_stop, index + 2., new_term, new_coeff, new_series_sum

    initial_term = a * torch.exp(
        -0.5 * h_squared * (1 - nega_squared)) / (2 * torch.pi)
    initial_sum = initial_term

    (_, _, _, _, series_sum) = tf.while_loop(
        cond=lambda stop, *_: tf.reduce_any(~stop),
        body=series_evaluation,
        loop_vars=(
            tf.cast(
                tf.zeros_like(h) * tf.zeros_like(a),
                dtype=tf.bool),
            tf.cast(3., dtype=dtype),
            initial_term,
            tf.ones_like(h) * tf.ones_like(a),
            initial_sum))
    return series_sum


def _owens_t_method5(h, a):
    """OwensT Method T5 which uses Gaussian Quadrature."""
    # Method T5, which is a gaussian quadrature approximation of the integral.

    # These are shifted and squared.
    quadrature_points = np.array([
        0.35082039676451715489E-02, 0.31279042338030753740E-01,
        0.85266826283219451090E-01, 0.16245071730812277011E+00,
        0.25851196049125434828E+00, 0.36807553840697533536E+00,
        0.48501092905604697475E+00, 0.60277514152618576821E+00,
        0.71477884217753226516E+00, 0.81475510988760098605E+00,
        0.89711029755948965867E+00, 0.95723808085944261843E+00,
        0.99178832974629703586E+00])
    quadrature_weights = np.array([
        0.18831438115323502887E-01, 0.18567086243977649478E-01,
        0.18042093461223385584E-01, 0.17263829606398753364E-01,
        0.16243219975989856730E-01, 0.14994592034116704829E-01,
        0.13535474469662088392E-01, 0.11886351605820165233E-01,
        0.10070377242777431897E-01, 0.81130545742299586629E-02,
        0.60419009528470238773E-02, 0.38862217010742057883E-02,
        0.16793031084546090448E-02])
    r = torch.square(a[..., tf.newaxis]) * quadrature_points
    log_integrand = -0.5 * torch.square(
        h[..., tf.newaxis]) * (1. + r) - torch.log1p(r)
    return torch.exp(torch.log(a) + torch.reduce_logsumexp(
        log_integrand + np.log(quadrature_weights), axis=-1))


def _owens_t_method6(h, a):
    # Method T6, which is a special case for when a is near 1.
    r = torch.atan2(1. - a, 1. + a)
    # When a = 1, T(h, 1) = 0.5 * ndtr(h) * (1 - ndtr(h)).
    # Thus, when a is close to 1, we add a correction term.
    normh = 0.5 * torch.erfc(h / np.sqrt(2.))
    result = 0.5 * normh * (1 - normh)
    return torch.where(
        torch.equal(r, 0.),
        result,
        result - r * torch.exp(
            -(1. - a) * torch.square(h) / (2 * r)) / (2 * torch.pi))


def _owens_t_regions(h, a):
    """Returns a list of Tensors describing the region of computation."""
    # We assume h >= 0, 0 <= a <= 1
    # Regions 1-7 that use T1.
    regions = []

    is_in_region1 = (h <= 0.06) & (a <= 0.025)
    is_in_region1 = is_in_region1 | (h <= 0.02) & (a <= 0.09)
    regions.append(is_in_region1)

    is_in_region2 = (h <= 0.02) & (a >= 0.09)
    is_in_region2 = (is_in_region2 |
                     (h >= 0.02) & (h <= 0.06) & (a >= 0.025) & (a <= 0.36))
    is_in_region2 = is_in_region2 | (h >= 0.06) & (h <= 0.09) & (a <= 0.09)
    regions.append(is_in_region2)

    is_in_region3 = (h >= 0.02) & (h <= 0.06) & (a >= 0.36)
    is_in_region3 = (is_in_region3 |
                     (h >= 0.06) & (h <= 0.09) & (a >= 0.09) & (a <= 0.5))
    is_in_region3 = (is_in_region3 |
                     (h >= 0.09) & (h <= 0.26) & (a >= 0.025) & (a <= 0.15))
    regions.append(is_in_region3)

    is_in_region4 = (h >= 0.06) & (h <= 0.125) & (a >= 0.9)
    regions.append(is_in_region4)

    is_in_region5 = (h >= 0.06) & (h <= 0.26) & (a >= 0.5) & (a <= 0.9)
    is_in_region5 = (is_in_region5 |
                     (h >= 0.09) & (h <= 0.26) & (a >= 0.15) & (a <= 0.5))
    is_in_region5 = (is_in_region5 |
                     (h >= 0.26) & (h <= 0.6) & (a >= 0.025) & (a <= 0.36))
    regions.append(is_in_region5)

    is_in_region6 = (h >= 0.26) & (h <= 0.6) & (a >= 0.36) & (a <= 0.9)
    is_in_region6 = is_in_region6 | (h >= 0.125) & (h <= 0.4) & (a >= 0.9)
    regions.append(is_in_region6)

    is_in_region7 = (h >= 0.6) & (h <= 1.7) & (a >= 0.15) & (a <= 0.36)
    regions.append(is_in_region7)

    is_in_region8 = (h >= 0.6) & (h <= 1.7) & (a >= 0.36) & (a <= 0.9)
    is_in_region8 = (is_in_region8 |
                     (h >= 0.4) & (h <= 1.6) & (a >= 0.9) & (a <= 0.99999))
    regions.append(is_in_region8)

    is_in_region9 = (h >= 4.8) & (a <= 0.09)
    regions.append(is_in_region9)

    is_in_region10 = (h >= 4.8) & (a >= 0.09) & (a <= 0.36)
    regions.append(is_in_region10)

    is_in_region11 = (h >= 4.8) & (a >= 0.36) & (a <= 0.5)
    regions.append(is_in_region11)

    is_in_region12 = (h >= 3.4) & (a >= 0.9)
    is_in_region12 = is_in_region12 | (h >= 3.36) & (a >= 0.36) & (a <= 0.9)
    is_in_region12 = is_in_region12 & ~is_in_region11
    regions.append(is_in_region12)

    is_in_region13 = (h >= 0.09) & (h <= 2.4) & (a <= 0.025)
    regions.append(is_in_region13)

    is_in_region14 = (h >= 0.6) & (h <= 1.7) & (a >= 0.025) & (a <= 0.09)
    regions.append(is_in_region14)

    is_in_region15 = (h >= 0.6) & (h <= 2.4) & (a >= 0.025) & (a <= 0.15)
    is_in_region15 = is_in_region15 & ~is_in_region14
    regions.append(is_in_region15)

    is_in_region16 = (h >= 1.7) & (h <= 2.4) & (a >= 0.15) & (a <= 0.36)
    is_in_region16 = is_in_region16 | (h >= 2.4) & (h <= 4.8) & (a <= 0.36)
    regions.append(is_in_region16)

    is_in_region17 = (h >= 1.6) & (h <= 3.4) & (a >= 0.9) & (a <= 0.99999)
    is_in_region17 = (is_in_region17 |
                      (h >= 1.7) & (h <= 3.4) & (a >= 0.36) & (a <= 0.9))
    regions.append(is_in_region17)

    # Near the line a = 1.
    is_in_region18 = (h >= 0.4) & (h <= 2.33) & (a >= 0.99999)
    regions.append(is_in_region18)

    return regions


def owens_t(h, a, name=None):
    # pylint: disable=line-too-long
    """Computes Owen's T function of `h` and `a` element-wise.

    Owen's T function is defined as the combined probability of the event `X > h`
    and `0 < Y < a * X`, where `X` and `Y` are independent standard normal
    random variables.

    In integral form this is defined as `1 / (2 * pi)` times the integral of
    `exp(-0.5 * h ** 2 * (1 + x ** 2)) / (1 + x ** 2)` from `0` to `a`.
    `h` and `a` can be any real number

    The Owen's T implementation below is based on
    ([Patefield and Tandy, 2000][1]).

    The Owen's T function has several notable properties which
    we list here for convenience. ([Owen, 1980][2], page 414)

    - P2.1  `T( h, 0)   =  0`
    - P2.2  `T( 0, a)   =  arctan(a) / (2 pi)`
    - P2.3  `T( h, 1)   =  Phi(h) (1 - Phi(h)) / 2`
    - P2.4  `T( h, inf) =  (1 - Phi(|h|)) / 2`
    - P2.5  `T(-h, a)   =  T(h, a)`
    - P2.6  `T( h,-a)   = -T(h, a)`
    - P2.7  `T( h, a) + T(a h, 1 / a) = Phi(h)/2 + Phi(ah)/2 - Phi(h) Phi(ah) - [a<0]/2`
    - P2.8  `T( h, a)   =  arctan(a)/(2 pi) - 1/(2 pi) int_0^h int_0^{ax}` exp(-(x**2 + y**2)/2) dy dx`
    - P2.9  `T( h, a)   =  arctan(a)/(2 pi) - int_0**h phi(x) Phi(a x) dx + Phi(h)/2 - 1/4`

    `[a<0]` uses Iverson bracket notation, i.e., `[a<0] = {1 if a<0 and 0 otherwise`.

    Let us also define P2.10 as:
    - P2.10  `T(inf, a) = 0`
    - Proof

      Note that result #10,010.6 ([Owen, 1980][2], pg 403) states that:
      `int_0^inf phi(x) Phi(a+bx) dx = Phi(a/rho)/2 + T(a/rho,b) where rho = sqrt(1+b**2).`
      Using `a=0`, this result is:
      `int_0^inf phi(x) Phi(bx) dx = 1/4 + T(0,b) = 1/4 + arctan(b) / (2 pi)`
      Combining this with P2.9 implies
      ```none
      T(inf, a)
       =  arctan(a)/(2 pi) - [ 1/4 + arctan(a) / (2 pi)]  + Phi(inf)/2 - 1/4
       = -1/4 + 1/2 -1/4 = 0.
      ```
      QED

    Args:
      h: A `float` `Tensor` defined as in `P({X > h, 0 < Y < a X})`. Must be
         broadcastable with `a`.
      a: A `float` `Tensor` defined as in `P({X > h, 0 < Y < a X})`. Must be
         broadcastable with `h`.
      name: A name for the operation (optional).

    Returns:
      owens_t: A `Tensor` with the same type as `h` and `a`,

    #### References

    [1]: Patefield, Mike, and D. A. V. I. D. Tandy. "Fast and accurate calculation
         of Owen’s T function." Journal of Statistical Software 5.5 (2000): 1-25.
         http://www.jstatsoft.org/v05/i05/paper
    [2]: Owen, Donald Bruce. "A table of normal integrals: A table."
         Communications in Statistics-Simulation and Computation 9.4 (1980):
         389-419.
    """
    # pylint: enable=line-too-long
    h = torch.tensor(h)
    a = torch.tensor(a)

    h = torch.abs(h)
    abs_a = torch.abs(a)
    # Remap arguments such that 0 <= a <= 1.
    modified_a = torch.where(
        abs_a <= 1.,
        abs_a,
        torch.reciprocal(abs_a))

    modified_h = torch.where(abs_a <= 1., h, abs_a * h)

    # For regions 1 - 8, we use method1 with different orders.

    regions = _owens_t_regions(modified_h, modified_a)

    # Short-circuit if we are not in the first 8 regions.
    order = 1.
    order = torch.where(regions[0], 2., order)
    order = torch.where(regions[1], 3., order)
    order = torch.where(regions[2], 4., order)
    order = torch.where(regions[3], 5., order)
    order = torch.where(regions[4], 7., order)
    order = torch.where(regions[5], 10., order)
    order = torch.where(regions[6], 12., order)
    order = torch.where(regions[7], 18., order)
    result = _owens_t_method1(modified_h, modified_a, order)

    # For regions 9, 10 and 11 we use method2 with different orders.
    order = 1.
    order = torch.where(regions[8], 10.), order)
    order = torch.where(regions[9], 20.), order)
    order = torch.where(regions[10], 30.), order)
    result = torch.where(
        regions[8] | regions[9] | regions[10],
        _owens_t_method2(modified_h, modified_a, order),
        result)

    # For region 12 we use method3.
    result = torch.where(
        regions[11], _owens_t_method3(modified_h, modified_a), result)

    # For regions 13, 14, 15 and 16 we use method4 with different orders.
    order = 1.)
    order = torch.where(regions[12], 4.), order)
    order = torch.where(regions[13], 7.), order)
    order = torch.where(regions[14], 8.), order)
    order = torch.where(regions[15], 20.), order)
    result = torch.where(
        regions[12] | regions[13] | regions[14] | regions[15],
        _owens_t_method4(modified_h, modified_a, order),
        result)

    # For region 17 we use method5.
    result = torch.where(
        regions[16], _owens_t_method5(modified_h, modified_a), result)

    # For region 18, we use method6.
    result = torch.where(
        regions[17], _owens_t_method6(modified_h, modified_a), result)

    result = torch.where(
        torch.equal(modified_h, 0.),
        torch.atan(modified_a) / (2 * torch.pi), result)

    # When a = 1, OwensT(h, 1) = ndtr(h) * (1 - ndtr(h))
    result = torch.where(
        torch.equal(modified_a, 1.),
        (0.125 * torch.erfc(-modified_h / math.sqrt(2.)) *
         torch.erfc(modified_h / math.sqrt(2.))), result)

    # When a = 0, we should return 0.
    result = torch.where(torch.equal(modified_a, 0.), 0.), result)

    normh = torch.erfc(h / math.sqrt(2.))
    normah = torch.erfc(abs_a * h / math.sqrt(2.))
    # Compensate for when |a| > 1.
    result = torch.where(
        abs_a > 1.,
        torch.where(
            abs_a * h <= 0.67,
            0.25 - 0.25 * torch.erf(
                h / math.sqrt(2.)) * torch.erf(abs_a * h / math.sqrt(2.)) - result,
            0.25 * (normh + normah - normh * normah) - result),
        result)

    result = torch.sign(a) * result

    result = torch.where(torch.isnan(a) | torch.isnan(h),
                         torch.nan,
                         result)
    return result
