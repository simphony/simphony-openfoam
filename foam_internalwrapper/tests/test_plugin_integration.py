import unittest


class TestPluginIntegration(unittest.TestCase):
    """Test case for Wrapper class"""

    def test_plugin_integration(self):
        """Test for plugin integration
        """

        # Assert that we can import the openfoam plugin
        from simphony.engine import openfoam_internal

        # Check that the expected top level objects are available
        self.assertTrue(hasattr(openfoam_internal, 'Wrapper'))


if __name__ == '__main__':
    unittest.main()
