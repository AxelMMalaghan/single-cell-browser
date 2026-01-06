from __future__ import annotations

__all__ = ["IDs", "reports_delete_id"]


class IDs:
    class Store:
        FILTER_STATE = "filter-state"
        USER_STATE = "user-state"
        SESSION_META = "session-metadata"
        ACTIVE_SESSION_ID = "active-session-id"
        ACTIVE_FIGURE_ID = "active-figure-id"

    class Control:
        # Core selectors
        DATASET_SELECT = "dataset-select"
        VIEW_SELECT = "view-select"

        CLUSTER_SELECT = "cluster-select"
        CONDITION_SELECT = "condition-select"
        SAMPLE_SELECT = "sample-select"
        CELLTYPE_SELECT = "celltype-select"
        GENE_SELECT = "gene-select"
        EMBEDDING_SELECT = "embedding-select"

        DIM_X = "dim-x-select"
        DIM_Y = "dim-y-select"
        DIM_Z_SELECT = "dim-z-select"  # The Dropdown
        DIM_Z_CONTAINER = "dim-z-container"

        OPTIONS_CHECKLIST = "options-checklist"
        COLOUR_SCALE_SELECT = "colour-scale-select"

        # Explore saved figures
        SAVED_FIGURE_SELECT = "saved-figure-select"
        SAVED_FIGURE_LOAD_BTN = "saved-figure-load-btn"
        FIGURE_LABEL_INPUT = "figure-label-input"
        SAVE_FIGURE_BTN = "save-figure-btn"
        SAVE_FIGURE_STATUS = "save-figure-status"

        # Explore sidebar / metadata
        SIDEBAR_DATASET_NAME = "sidebar-dataset-name"
        SIDEBAR_DATASET_META = "sidebar-dataset-meta"

        # Explore filters (containers)
        CLUSTER_FILTER_CONTAINER = "cluster-filter-container"
        CONDITION_FILTER_CONTAINER = "condition-filter-container"
        SAMPLE_FILTER_CONTAINER = "sample-filter-container"
        CELLTYPE_FILTER_CONTAINER = "celltype-filter-container"
        GENE_FILTER_CONTAINER = "gene-filter-container"
        EMBEDDING_FILTER_CONTAINER = "embedding-filter-container"
        DIM_FILTER_CONTAINER = "dim-filter-container"
        OPTIONS_CONTAINER = "options-container"

        # Graph + downloads
        MAIN_GRAPH = "main-graph"
        DOWNLOAD_DATA = "download-data"
        DOWNLOAD_DATA_BTN = "download-data-btn"

        # Status bar
        STATUS_BAR = "status-bar"

        # Dataset manager / preview
        DM_OBS_PREVIEW = "dm-obs-preview"
        DM_STATUS_TEXT = "dm-status-text"
        DM_SUMMARY_TEXT = "dm-summary-text"
        DM_CLUSTER_KEY = "dm-cluster-key"
        DM_CONDITION_KEY = "dm-condition-key"
        DM_SAMPLE_KEY = "dm-sample-key"
        DM_CELLTYPE_KEY = "dm-celltype-key"
        DM_EMBEDDING_KEY = "dm-embedding-key"
        DM_CURRENT_DATASET = "dm-current-dataset"

        DM_SAVE_STATUS = "dm-save-status"
        DM_SAVE_BTN = "dm-save-btn"

        DM_IMPORT_STATUS = "dm-import-status"
        DM_UPLOAD = "dm-upload"

        # Reports
        REPORTS_EXPORT_BTN = "reports-metadata_io-btn"
        REPORTS_UPLOAD = "reports-upload"
        REPORTS_SUMMARY_TEXT = "reports-summary-text"
        REPORTS_FIGURE_LIST = "reports-figure-list"
        REPORTS_DOWNLOAD_SESSION = "reports-download-session"
        REPORTS_IMPORT_BANNER = "reports-import-banner"

    class Pattern:
        # pattern-matching "type" strings
        REPORTS_DELETE = "reports-delete-figure"


def reports_delete_id(fig_id: str) -> dict:
    return {"type": IDs.Pattern.REPORTS_DELETE, "index": fig_id}
