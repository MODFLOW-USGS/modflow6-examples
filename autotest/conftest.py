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
def run(request) -> bool:
    return request.config.getoption("--run")


@pytest.fixture(scope="session")
def plot(request) -> bool:
    return request.config.getoption("--plot")


@pytest.fixture(scope="session")
def save(request) -> bool:
    return request.config.getoption("--save")


def pytest_addoption(parser):
    parser.addoption("--run", action="store_true", default=False, help="Whether to run models")
    parser.addoption("--plot", action="store_true", default=False, help="Whether to plot results (requires --run)")
    parser.addoption("--save", action="store_true", default=False, help="Whether to save plots (requires --run and --plot)")


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
