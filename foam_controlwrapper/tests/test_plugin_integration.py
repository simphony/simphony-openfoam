import unittest


class TestPluginIntegration(unittest.TestCase):
    """Test case for Wrapper class"""

    def test_plugin_integration(self):
        """Test for plugin integration
        """

        # Assert that we can import the openfoam plugin
        from simphony.engine import openfoam_file_io

        # Check that the expected top level objects are available
        self.assertTrue(hasattr(openfoam_file_io, 'Wrapper'))


if __name__ == '__main__':
    unittest.main()
