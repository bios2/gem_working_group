import numpy as np

from gem.vegetation import logistic_growth_delta


def test_logistic_growth_zero_at_carrying_capacity():
    B = np.array([100.0, 100.0, 100.0])
    r = np.array([0.1, 0.2, 0.3])
    K = np.array([100.0, 100.0, 100.0])
    delta = logistic_growth_delta(B, r, K, dt=1.0)
    np.testing.assert_allclose(delta, 0.0)


def test_logistic_growth_runs_on_grid():
    shape = (4, 5, 3)
    B = np.full(shape, 50.0)
    r = np.full(shape, 0.1)
    K = np.full(shape, 100.0)
    delta = logistic_growth_delta(B, r, K, dt=1.0)
    assert delta.shape == shape
    assert np.all(delta > 0)
