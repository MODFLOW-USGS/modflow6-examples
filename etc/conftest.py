from pathlib import Path
import pytest

PROJ_ROOT = Path(__file__).parent.parent

@pytest.fixture(scope="session")
def run(pytestconfig):
    return pytestconfig.getoption("run")


def pytest_addoption(parser):
    parser.addoption("--run", action="store", default=False)


def pytest_generate_tests(metafunc):
    if "example" in metafunc.fixturenames:
        run = metafunc.config.getoption("run")
        scripts = {
            file.name: (run, file)
            for file in sorted((PROJ_ROOT / "scripts").glob("ex-*.py"))
        }
        metafunc.parametrize(
            "example", scripts.values(), ids=scripts.keys()
        )
