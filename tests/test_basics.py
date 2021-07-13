import pytest
from subprocess import check_call

@pytest.mark.parametrize("command", ["-version", "-this-is-not-a-valid-parameter"] )
def test_basis(command):
    if command == "-this-is-not-a-valid-parameter":
        with pytest.raises(Exception):
            check_call(["wsclean", command])
    else:
        check_call(["wsclean", command])