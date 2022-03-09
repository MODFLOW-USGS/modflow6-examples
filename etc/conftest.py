import os


def pytest_addoption(parser):
    parser.addoption("--run", action="store", default=False)


def pytest_generate_tests(metafunc):
    if "build" in metafunc.fixturenames:
        run = metafunc.config.getoption("run")
        scripts_dict = {
            file_name: (run, file_name)
            for file_name in sorted(os.listdir(os.path.join("..", "scripts")))
            if file_name.endswith(".py") and file_name.startswith("ex-")
        }
        metafunc.parametrize(
            "build", scripts_dict.values(), ids=scripts_dict.keys()
        )
