import unittest
import Newton_Raphson as nr


class TestNewtonRaphson(unittest.TestCase):

    def test_good_enough(self):
        self.assertFalse(nr.good_enough(2, 1))
        self.assertTrue(nr.good_enough(1, 1.00000001))


if __name__ == '__main__':
    unittest.main()
