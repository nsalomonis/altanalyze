import unittest

from genapi import GenCloud

class TestLogin(unittest.TestCase):

    def test_login(self):
        GenCloud('admin', 'admin', 'http://gendev:10180')

if __name__ == '__main__':
    unittest.main()
