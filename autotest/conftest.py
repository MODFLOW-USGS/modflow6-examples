from pathlib import Path
import pytest

PROJ_ROOT = Path(__file__).parents[1]
SCRIPTS_PATH = PROJ_ROOT / "scripts"
IMAGES_PATH = PROJ_ROOT / "images"
FIGURES_PATH = PROJ_ROOT / "figures"
EXAMPLES_PATH = PROJ_ROOT / "examples"
NOTEBOOKS_PATH = PROJ_ROOT / "notebooks"
RTD_PATH = PROJ_ROOT / ".doc" / "_notebooks"


@pytest.fixture(scope="session")
def run(pytestconfig):
    return pytestconfig.getoption("run")


def pytest_addoption(parser):
    parser.addoption("--run", action="store", default=False)


def pytest_generate_tests(metafunc):
    # make directories if needed
    IMAGES_PATH.mkdir(exist_ok=True)
    FIGURES_PATH.mkdir(exist_ok=True)
    EXAMPLES_PATH.mkdir(exist_ok=True)
    NOTEBOOKS_PATH.mkdir(exist_ok=True)
    RTD_PATH.mkdir(exist_ok=True, parents=True)

    # generate example scenarios
    if "example" in metafunc.fixturenames:
        run = metafunc.config.getoption("run")
        scripts = {
            file.name: (run, file)
            for file in sorted(SCRIPTS_PATH.glob("ex-*.py"))
        }
        metafunc.parametrize(
            "example", scripts.values(), ids=scripts.keys()
        )
