from .. import stats


def test_concentrated_data_produces_entropy_0():
    entropy = stats.base4_entropy([1, 0, 0, 0], axis=0)
    assert entropy == 0


def test_concentrated_unnormalized_data_produces_entropy_0():
    entropy = stats.base4_entropy([1000, 0, 0, 0], axis=0)
    assert entropy == 0


def test_balanced_data_produces_entropy_1():
    entropy = stats.base4_entropy([0.25, 0.25, 0.25, 0.25], axis=0)
    assert entropy == 1


def test_balanced_unnormalized_data_produces_entropy_1():
    entropy = stats.base4_entropy([20, 20, 20, 20], axis=0)
    assert entropy == 1
