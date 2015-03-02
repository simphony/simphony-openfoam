import unittest


class TestPluginIntegration(unittest.TestCase):
    """Test case for FoamControlWrapper class"""

    def test_plugin_integration(self):
        """Test to run pitzDaily example in OpenFoam"
        """

        # Assert that we can import the openfoam plugin
        from simphony.engine import openfoam

        # Check that the expected top level objects are available
        self.assertTrue(hasattr(openfoam, 'FoamControlWrapper'))


if __name__ == '__main__':
    unittest.main()
