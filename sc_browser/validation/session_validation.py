from __future__ import annotations

from typing import Any

from sc_browser.validation.errors import ValidationIssue, ValidationError


def validate_session_import_dict(obj: Any) -> None:
    """
    Validate the raw JSON dict BEFORE you build a SessionMetadata.
    This prevents half-valid imports from poisoning app state.
    """
    issues: list[ValidationIssue] = []

    if not isinstance(obj, dict):
        raise ValidationError([ValidationIssue("SESSION_TYPE", "Uploaded session metadata must be a JSON object.")])

    figs = obj.get("figures", [])
    if figs is None:
        figs = []

    if not isinstance(figs, list):
        issues.append(ValidationIssue("SESSION_FIGURES_TYPE", "figures must be a list."))

    # Minimal checks per figure (cheap + effective)
    if isinstance(figs, list):
        for i, f in enumerate(figs):
            if not isinstance(f, dict):
                issues.append(ValidationIssue("FIG_TYPE", f"figures[{i}] must be an object."))
                continue
            if not f.get("dataset_key"):
                issues.append(ValidationIssue("FIG_DATASET_KEY", f"figures[{i}].dataset_key missing."))
            if not f.get("view_id"):
                issues.append(ValidationIssue("FIG_VIEW_ID", f"figures[{i}].view_id missing."))
            fs = f.get("filter_state") or f.get("state") or {}
            if fs and not isinstance(fs, dict):
                issues.append(ValidationIssue("FIG_FILTER_STATE", f"figures[{i}].filter_state/state must be an object."))

    if issues:
        raise ValidationError(issues)
