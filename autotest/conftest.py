from pathlib import Path

import pytest

PROJ_ROOT = Path(__file__).parents[1]
SCRIPTS_PATH = PROJ_ROOT / "scripts"
TABLES_PATH = PROJ_ROOT / "tables"
IMAGES_PATH = PROJ_ROOT / "images"
FIGURES_PATH = PROJ_ROOT / "figures"
EXAMPLES_PATH = PROJ_ROOT / "examples"
NOTEBOOKS_PATH = PROJ_ROOT / "notebooks"
RTD_PATH = PROJ_ROOT / ".doc" / "_notebooks"
EXCLUDE = []


@pytest.fixture(scope="session")
def write(request) -> bool:
    return not request.config.getoption("--no-write")


@pytest.fixture(scope="session")
def run(request) -> bool:
    return not request.config.getoption("--no-run")


@pytest.fixture(scope="session")
def plot(request) -> bool:
    return not request.config.getoption("--no-plot")


@pytest.fixture(scope="session")
def gif(request) -> bool:
    return not request.config.getoption("--no-gif")


def pytest_addoption(parser):
    parser.addoption("--no-write", action="store_true", default=False, help="Disable model input file creation")
    parser.addoption("--no-run", action="store_true", default=False, help="Disable model runs")
    parser.addoption("--no-plot", action="store_true", default=False, help="Disable static plot creation")
    parser.addoption("--no-gif", action="store_true", default=False, help="Disable GIF creation")


def pytest_generate_tests(metafunc):
    # make directories if needed
    TABLES_PATH.mkdir(exist_ok=True)
    IMAGES_PATH.mkdir(exist_ok=True)
    FIGURES_PATH.mkdir(exist_ok=True)
    EXAMPLES_PATH.mkdir(exist_ok=True)
    NOTEBOOKS_PATH.mkdir(exist_ok=True)
    RTD_PATH.mkdir(exist_ok=True, parents=True)

    # generate example scenarios
    if "example_script" in metafunc.fixturenames:
        scripts = {
            file.name: file
            for file in sorted(SCRIPTS_PATH.glob("ex-*.py"))
            if file.stem not in EXCLUDE
        }
        metafunc.parametrize("example_script", scripts.values(), ids=scripts.keys())
