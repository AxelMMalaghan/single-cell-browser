# UI Structure and Styling Reference

This document captures the current UI structure and styling for the single-cell browser Dash app.
It is intended for future agentic use when navigating or modifying the UI.

## UI Structure

- App layout is built in `sc_browser/ui/layout/build_layout.py` with a sticky navbar, app-level
  `dcc.Store`s, and a `dcc.Tabs` scaffold for Dataset Importer, Dataset Preview, Explore, and Reports.
  Explore uses a two-column grid (md=3 sidebar for View/Label + Filters, md=9 plot panel).
- Navbar includes logo + title/subtitle and a dataset selector block with active dataset info; all
  layout/IDs live in `sc_browser/ui/layout/build_navbar.py`.
- Sidebar cards are split between “View & label” (saved view dropdown, view type, figure label,
  load/save buttons) in `sc_browser/ui/layout/build_view_and_label_panel.py` and “Filters” (dataset
  summary + filter dropdowns + options checklist) in `sc_browser/ui/layout/build_filter_panel.py`.
- Main plot panel is a card with a graph + download action in
  `sc_browser/ui/layout/build_plot_panel.py`.
- Dataset Importer is a two-column card layout: current dataset status + upload on the left, mapping
  dropdowns and save on the right in `sc_browser/ui/layout/build_dataset_import_panel.py`.
- Dataset Preview is a full-width card using the maincard styles in
  `sc_browser/ui/layout/build_dataset_preview_panel.py`.
- Reports is a card with summary + export/import actions and a figure list table in
  `sc_browser/ui/layout/build_reports_panel.py`.

## Styling

- Global look is Bootstrap “FLATLY” plus `sc_browser/ui/assets/style.css` (loaded via
  `assets_folder` in `sc_browser/ui/dash_app.py`), giving a light gray canvas with strong, dark
  navbar contrast.
- Typography: system-ui stack, uppercase/letter-spaced nav and tab labels, small UI text sizes for
  form labels and tables.
- Navbar: full-width sticky banner with deep navy background, heavy shadow, and dataset selector
  block aligned right.
- Cards: borderless by default, but custom cards (`.scb-sidebar`, `.scb-viewlabel-card`,
  `.scb-maincard`) add light borders, soft shadows, and subtle header backgrounds.
- Controls: dropdowns and inputs are flattened with small font sizes, light borders, and explicit
  z-index fixes for dropdown overflow.
- Buttons: primary uses dark theme with pronounced shadow; secondary is light with subtle
  border/shadow.
- Tabs: uppercase, minimal chrome; selected tab gets a white background and bottom border.
- Utility styles for uploads, reports tables, and Dash DataTable font consistency are centralized
  in `sc_browser/ui/assets/style.css`.

## Layout Parameters (Exact Values)

Values below are taken from `sc_browser/ui/assets/style.css` and layout files.

### Global and Containers

- `body` background: `#f3f4f6`, text: `#111827`.
- Root container `.scb-root`: `padding: 0 1.5rem 2rem`, no max-width, no margin.
- Bootstrap cards default border removed: `.card { border: none; }`.

### Navbar

- `.scb-navbar`:
  - `padding: 1.25rem 1.5rem`
  - `margin: 0 -1.5rem 0` (full-bleed over root padding)
  - `min-height: 130px`
  - `background-color: #020617`
  - `position: sticky; top: 0; z-index: 100`
  - `box-shadow: 0 8px 24px rgba(0,0,0,0.4), 0 0 0 1px rgba(0,0,0,0.7)`
- Title `h2` in navbar:
  - `font-size: 1.7rem`
  - `font-weight: 720`
  - `letter-spacing: 0.08em`
  - `text-transform: uppercase`
- Subtitle `#navbar-subtitle`:
  - `font-size: 0.85rem`
  - `opacity: 0.9`
- Dataset block `.navbar-dataset-block`:
  - `padding: 0.5rem 0.75rem`
  - `border-radius: 0.5rem`
  - `minWidth: 280px; maxWidth: 380px; marginRight: 24px` (inline style in `build_navbar.py`)
- Dataset title `.navbar-dataset-title`:
  - `font-size: 0.8rem`, `letter-spacing: 0.12em`, `text-align: right`
- Dataset subtitle `.navbar-dataset-subtitle`:
  - `font-size: 0.75rem`, `text-align: right`
- Navbar dataset dropdown `.scb-navbar .scb-dataset-dropdown .Select-control`:
  - `font-size: 0.85rem`
  - `min-height: 36px`
  - `height: 36px`
  - Text vertical alignment uses `line-height: 34px` on value/placeholder.

### Tabs

- Tabs container `#page-tabs`: `margin-top: 1rem`.
- Tab `.tab`:
  - `padding: 0.45rem 1rem`
  - `margin-right: 0.3rem`
  - `border-bottom: 2px solid transparent`
  - `font-size: 0.8rem`
  - `letter-spacing: 0.08em`
- Selected tab `.tab--selected`:
  - `border-bottom: 2px solid #111827`
  - `box-shadow: 0 -2px 10px rgba(15, 23, 42, 0.08)`

### Cards and Panels

- Sidebar cards `.scb-sidebar`, `.scb-viewlabel-card`:
  - `border-radius: 4px`
  - `border: 1px solid #e5e7eb`
  - `box-shadow: 0 8px 18px rgba(15, 23, 42, 0.06)`
- Sidebar/view header `.card-header` inside sidebar/view:
  - `padding: 0.55rem 0.9rem`
  - `font-size: 0.8rem`
  - `letter-spacing: 0.08em`
- Sidebar/view body `.card-body`:
  - `padding: 0.85rem 0.9rem 0.9rem`
- Main plot card `.scb-maincard`:
  - `border-radius: 4px`
  - `border: 1px solid #e5e7eb`
  - `box-shadow: 0 12px 28px rgba(15, 23, 42, 0.08)`
- Main card header `.scb-maincard .card-header`:
  - `padding: 0.55rem 0.9rem`
- Main card body `.scb-main-body`:
  - `padding: 0.75rem 0.85rem 0.9rem`
- Graph `#main-graph`:
  - `border-radius: 4px`
  - `border: 1px solid #e5e7eb`
  - `background: #ffffff`
  - inline height in `build_plot_panel.py`: `style={"height": "650px"}`

### Controls and Text

- Form labels `.form-label`:
  - `font-size: 0.78rem`
  - `font-weight: 600`
  - `margin-bottom: 0.25rem`
- Text inputs `input.form-control`:
  - `border-radius: 4px`
  - `font-size: 0.8rem`
  - `border: 1px solid #d1d5db`
- Dropdowns `.Select-control` within `.scb-root`:
  - `border-radius: 4px`
  - `min-height: 2.1rem`
  - `font-size: 0.8rem`
- Small text `.small, .small *`: `font-size: 0.78rem`.

### Buttons

- Primary `.btn-primary`:
  - `font-size: 0.8rem`
  - `padding: 0.25rem 0.9rem`
  - `border-radius: 4px`
  - `box-shadow: 0 8px 18px rgba(15, 23, 42, 0.3)`
- Secondary `.btn-secondary`:
  - `font-size: 0.8rem`
  - `padding: 0.25rem 0.85rem`
  - `border-radius: 4px`
  - `box-shadow: 0 4px 10px rgba(15, 23, 42, 0.08)`
