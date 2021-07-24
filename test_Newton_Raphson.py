import unittest
import Newton_Raphson as nr


class TestNewtonRaphson(unittest.TestCase):

    def test_good_enough(self):
        self.assertFalse(nr.good_enough(2, 1))


if __name__ == '__main__':
    unittest.main()
