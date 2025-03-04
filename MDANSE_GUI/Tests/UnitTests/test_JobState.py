import pytest

from MDANSE_GUI.Tabs.Models.JobHolder import JobEntry


@pytest.fixture(scope="module")
def temporary_jobentry() -> JobEntry:
    return JobEntry()

@pytest.mark.parametrize("task, result", [
    ("start", "Running"),
    ("fail", "Failed"),
])
def test_task(temporary_jobentry: JobEntry, task, result):
    getattr(temporary_jobentry._current_state, task)()
    assert temporary_jobentry._current_state._label == result
