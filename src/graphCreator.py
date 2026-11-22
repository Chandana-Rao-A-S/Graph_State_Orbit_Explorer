"""Graph State Orbit Explorer: an interactive editor for graphs and their LC orbits.

Run from the project root with::

    python -m src.graphCreator

The main window is a graph editor with a live statistics panel. Local
complementation (LC) can be applied at any vertex, the LC orbit of the current
graph can be computed and browsed, and the orbit can be drawn as a graph of
graphs in which two members are joined when one LC maps one to the other.

This module contains only the Tk interface. The mathematics lives in
``graphAnalytics`` (invariants, Schmidt ranks), ``graphOperations`` (LC and
orbit search) and ``toolset`` (canonical labelling).
"""

import csv
import io
import math
import multiprocessing
import threading
import tkinter as tk
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
from tkinter import filedialog, messagebox, simpledialog, ttk

import networkx as nx

from src import graphOperations
from src.graphAnalytics import (
    ORBIT_TABLE_COLUMNS,
    graph_stats,
    greedy_coloring,
    is_k_colorable,
    orbit_table_row,
    schmidt_rank_profile,
)
from src.toolset import canonical_relabel_n as canonical_relabel


# Largest graph for which the stats panel enumerates all bipartitions on every edit.
LIVE_SCHMIDT_MAX_NODES = 16
# Largest graph accepted by the "Compute Schmidt" button.
BUTTON_SCHMIDT_MAX_NODES = 20
# Orbits with more members than this are drawn only as neighbourhoods.
OVERVIEW_FULL_MAX_NODES = 500
# Orbits with at least this many members get their statistics table in parallel.
PARALLEL_MIN_ROWS = 200

BG = "#2b2b2b"
PALETTE = ["#4e79a7", "#f28e2b", "#e15759", "#76b7b2", "#59a14f",
           "#edc949", "#af7aa1", "#ff9da7", "#9c755f", "#bab0ab"]


# ─────────────────────────────────────────────────────────────────────────────
#  Help texts
# ─────────────────────────────────────────────────────────────────────────────

STATS_HELP = {
    "nodes": "Number of vertices n (qubits). LC invariant.",
    "edges": "Number of edges m.",
    "density": "Edge density m / (n choose 2), between 0 and 1.",
    "clustering": "Average clustering coefficient (Watts–Strogatz).",
    "assortativity": "Degree assortativity: Pearson correlation of the degrees at the two "
                     "ends of an edge. Undefined (—) for regular graphs.",
    "connected": "Whether the graph is connected. LC invariant: LC never joins or splits "
                 "connected components.",
    "components": "Number of connected components. LC invariant.",
    "diameter_giant": "Diameter of the largest connected component.",
    "avg_path_len_giant": "Average shortest-path length within the largest connected component.",
    "bipartite": "Whether the graph is 2-colourable (no odd cycles). Not an LC invariant.",
    "bipartition_sizes": "Sizes of the two colour classes of a connected bipartite graph.",
    "planar": "Whether the graph can be drawn in the plane without crossing edges.",
    "clique_number": "ω(G): size of a largest complete subgraph.",
    "chromatic_lower_bound": "Lower bound on the chromatic number χ(G): ω(G) ≤ χ(G).",
    "chromatic_upper_bound": "Upper bound on χ(G): colours used by the greedy DSATUR colouring.",
    "deg_min": "Minimum vertex degree.",
    "deg_avg": "Average vertex degree.",
    "deg_max": "Maximum vertex degree.",
    "triangles": "Number of triangles (3-cycles).",
    "schmidt_rank_max": "Largest Schmidt rank of the graph state |G⟩ over all bipartitions A|B. "
                        "Across A|B the Schmidt rank is 2^r, with r the GF(2) rank of the "
                        "adjacency matrix between A and B. LC invariant. Computed live for "
                        f"n ≤ {LIVE_SCHMIDT_MAX_NODES}.",
    "_schmidt_tv": "Number of bipartitions A|B with each Schmidt rank. All 2^(n−1) − 1 "
                   "bipartitions are counted. The whole distribution is an LC invariant.",
    "_orbit_size": "Number of pairwise non-isomorphic graphs in the LC orbit.",
}

ORBIT_STATS_COL_HELP = {
    "graph6": "Canonical graph6 string. Two graphs are isomorphic if and only if their "
              "strings are equal.",
    "n": "Number of vertices n. LC invariant.",
    "m": "Number of edges m.",
    "connected": "Whether the graph is connected. LC invariant.",
    "diam(Gc)": "Diameter of the largest connected component.",
    "planar": "Whether the graph is planar.",
    "bipartite": "Whether the graph is bipartite. Not an LC invariant.",
    "omega": "Clique number ω(G).",
    "chi_greedy": "Colours used by the greedy DSATUR colouring (upper bound on χ).",
    "Schmidt_max": "Largest Schmidt rank over all bipartitions. LC invariant, so equal "
                   "along the whole orbit.",
}


# ─────────────────────────────────────────────────────────────────────────────
#  Tooltips
# ─────────────────────────────────────────────────────────────────────────────

class Tooltip:
    """Dark tooltip that appears after the pointer rests on a widget."""

    def __init__(self, widget, text: str, delay: int = 400):
        self.widget = widget
        self.text = text
        self.delay = delay
        self._after_id = None
        self._tip = None
        widget.bind("<Enter>", self._schedule, add="+")
        widget.bind("<Leave>", self._hide, add="+")
        widget.bind("<Motion>", self._move, add="+")

    def _schedule(self, _event):
        self._cancel()
        self._after_id = self.widget.after(self.delay, self._show)

    def _cancel(self):
        if self._after_id:
            self.widget.after_cancel(self._after_id)
            self._after_id = None

    def _show(self):
        if self._tip or not self.text:
            return
        x = self.widget.winfo_rootx() + 16
        y = self.widget.winfo_rooty() + self.widget.winfo_height() + 4
        self._tip = _tooltip_window(self.widget, self.text, x, y)

    def _move(self, event):
        if self._tip:
            self._tip.wm_geometry(f"+{event.x_root + 12}+{event.y_root + 12}")

    def _hide(self, _event=None):
        self._cancel()
        if self._tip:
            self._tip.destroy()
            self._tip = None


def _tooltip_window(parent, text, x, y):
    tip = tk.Toplevel(parent)
    tip.wm_overrideredirect(True)
    tip.attributes("-topmost", True)
    frame = tk.Frame(tip, background="#333", borderwidth=1, relief="solid")
    frame.pack()
    tk.Label(frame, text=text, justify="left", background="#333", foreground="#fff",
             padx=8, pady=6, wraplength=360).pack()
    tip.wm_geometry(f"+{x}+{y}")
    return tip


def add_treeview_header_tooltips(tree: ttk.Treeview, help_map: dict, delay_ms: int = 350):
    """Show a tooltip while the pointer rests on a Treeview column header."""
    state = {"tip": None, "timer": None}

    def hide():
        if state["timer"]:
            tree.after_cancel(state["timer"])
            state["timer"] = None
        if state["tip"]:
            state["tip"].destroy()
            state["tip"] = None

    def show(text, x, y):
        hide()
        state["tip"] = _tooltip_window(tree, text, x, y)

    def on_motion(event):
        if tree.identify_region(event.x, event.y) != "heading":
            hide()
            return
        column_id = tree.identify_column(event.x)
        if not column_id[1:].isdigit():
            hide()
            return
        index = int(column_id[1:]) - 1
        columns = tree["columns"]
        text = help_map.get(columns[index], "") if 0 <= index < len(columns) else ""
        hide()
        if text:
            state["timer"] = tree.after(
                delay_ms, lambda: show(text, event.x_root + 12, event.y_root + 12)
            )

    tree.bind("<Motion>", on_motion, add="+")
    tree.bind("<Leave>", lambda _e: hide(), add="+")


# ─────────────────────────────────────────────────────────────────────────────
#  Drawing primitives
# ─────────────────────────────────────────────────────────────────────────────

class Node:
    """A vertex drawn on a GraphView canvas."""

    RADIUS = 24
    FILL = "#007acc"
    OUTLINE = "black"
    HIGHLIGHT = "#ffb347"
    MARK = "#00e676"
    FONT = ("Helvetica", 11, "bold")

    def __init__(self, view, x: float, y: float, label):
        self.view = view
        self.x = x
        self.y = y
        self.label = str(label)
        self.edges = set()
        self.selected = False
        self.marked = False
        self._hovered = False

        r = self.RADIUS
        canvas = view.canvas
        self.oval_id = canvas.create_oval(x - r, y - r, x + r, y + r, fill=self.FILL,
                                          outline=self.OUTLINE, width=3, tags=("node",))
        self.text_id = canvas.create_text(x, y, text=self.label, fill="white",
                                          font=self.FONT, tags=("node",))
        for item in (self.oval_id, self.text_id):
            view.item_to_node[item] = self
            canvas.tag_bind(item, "<Enter>", lambda _e: self._set_hover(True))
            canvas.tag_bind(item, "<Leave>", lambda _e: self._set_hover(False))
        view.nodes.append(self)

    def _set_hover(self, on: bool):
        self._hovered = on
        self._refresh_outline()

    def _refresh_outline(self):
        if self.marked:
            outline, width = self.MARK, 6
        elif self.selected or self._hovered:
            outline, width = self.HIGHLIGHT, 3
        else:
            outline, width = self.OUTLINE, 3
        self.view.canvas.itemconfig(self.oval_id, outline=outline, width=width)

    def set_selected(self, on: bool):
        self.selected = on
        self._refresh_outline()

    def set_marked(self, on: bool):
        self.marked = on
        self._refresh_outline()
        if on:
            self.view.canvas.tag_raise(self.oval_id)
            self.view.canvas.tag_raise(self.text_id)

    def set_fill(self, color: str):
        self.view.canvas.itemconfig(self.oval_id, fill=color)

    def set_label(self, label: str):
        self.label = label
        self.view.canvas.itemconfig(self.text_id, text=label)

    def move_to(self, x: float, y: float):
        dx, dy = x - self.x, y - self.y
        self.x, self.y = x, y
        self.view.canvas.move(self.oval_id, dx, dy)
        self.view.canvas.move(self.text_id, dx, dy)
        for edge in self.edges:
            edge.update_position()

    def delete(self):
        for edge in list(self.edges):
            edge.delete()
        self.view.canvas.delete(self.oval_id)
        self.view.canvas.delete(self.text_id)
        self.view.item_to_node.pop(self.oval_id, None)
        self.view.item_to_node.pop(self.text_id, None)
        self.view.nodes.remove(self)


class Edge:
    """An undirected edge drawn on a GraphView canvas."""

    COLOR = "#f5f5f5"
    HIGHLIGHT = "#ffd166"

    def __init__(self, view, node_a: Node, node_b: Node):
        if node_a is node_b:
            raise ValueError("Self-loops are not allowed in graph states.")
        if view.edge_between(node_a, node_b) is not None:
            raise ValueError("These vertices are already joined by an edge.")
        self.view = view
        self.node_a = node_a
        self.node_b = node_b
        self.selected = False
        self.line_id = view.canvas.create_line(*self._coords(), width=3, fill=self.COLOR,
                                               tags=("edge",))
        view.canvas.tag_lower(self.line_id)
        view.canvas.tag_bind(self.line_id, "<Enter>", lambda _e: self._style(True))
        view.canvas.tag_bind(self.line_id, "<Leave>", lambda _e: self._style(False))
        node_a.edges.add(self)
        node_b.edges.add(self)
        view.item_to_edge[self.line_id] = self

    def other(self, node: Node) -> Node:
        return self.node_b if node is self.node_a else self.node_a

    def _coords(self):
        """End points on the borders of the two vertex circles."""
        dx = self.node_b.x - self.node_a.x
        dy = self.node_b.y - self.node_a.y
        dist = math.hypot(dx, dy) or 1.0
        ux, uy = dx / dist, dy / dist
        r = Node.RADIUS
        return (self.node_a.x + ux * r, self.node_a.y + uy * r,
                self.node_b.x - ux * r, self.node_b.y - uy * r)

    def _style(self, hover: bool):
        on = hover or self.selected
        self.view.canvas.itemconfig(self.line_id, width=5 if on else 3,
                                    fill=self.HIGHLIGHT if on else self.COLOR)

    def set_selected(self, on: bool):
        self.selected = on
        self._style(False)

    def update_position(self):
        self.view.canvas.coords(self.line_id, *self._coords())

    def delete(self):
        self.view.canvas.delete(self.line_id)
        self.node_a.edges.discard(self)
        self.node_b.edges.discard(self)
        self.view.item_to_edge.pop(self.line_id, None)


class GraphView(ttk.Frame):
    """A canvas that holds a drawn graph and converts it to and from NetworkX."""

    def __init__(self, parent):
        super().__init__(parent)
        self.canvas = tk.Canvas(self, bg=BG, highlightthickness=0)
        self.canvas.pack(fill=tk.BOTH, expand=True)
        self.nodes: list[Node] = []
        self.item_to_node: dict[int, Node] = {}
        self.item_to_edge: dict[int, Edge] = {}
        self._size = None
        self.canvas.bind("<Configure>", self._on_resize)

    # ----- queries -----
    def edges(self) -> list[Edge]:
        return list(self.item_to_edge.values())

    def edge_between(self, a: Node, b: Node):
        return next((e for e in a.edges if e.other(a) is b), None)

    def item_at(self, x: float, y: float, halo: int = 4):
        """The vertex at (x, y) if there is one, otherwise an edge near (x, y), or None."""
        hits = self.canvas.find_overlapping(x - halo, y - halo, x + halo, y + halo)
        for item in reversed(hits):
            if item in self.item_to_node:
                return self.item_to_node[item]
        for item in reversed(hits):
            if item in self.item_to_edge:
                return self.item_to_edge[item]
        return None

    def to_networkx(self, use_labels: bool = False) -> nx.Graph:
        """The drawn graph. Vertices are 0..n-1 in creation order, or the labels."""
        index = {node: (node.label if use_labels else i) for i, node in enumerate(self.nodes)}
        G = nx.Graph()
        G.add_nodes_from(index[node] for node in self.nodes)
        G.add_edges_from((index[e.node_a], index[e.node_b]) for e in self.edges())
        return G

    # ----- drawing -----
    def clear(self):
        for node in list(self.nodes):
            node.delete()

    def _canvas_size(self):
        self.canvas.update_idletasks()
        w, h = self.canvas.winfo_width(), self.canvas.winfo_height()
        if w <= 1 or h <= 1:
            w, h = self.canvas.winfo_reqwidth(), self.canvas.winfo_reqheight()
        return w, h

    def draw(self, G: nx.Graph, labels: dict | None = None, scale: float = 0.8) -> dict:
        """Replace the drawing by ``G`` with its vertices on a circle.

        Returns a mapping from the vertices of ``G`` to the drawn Node objects.
        """
        self.clear()
        vertices = list(G.nodes)
        n = len(vertices)
        w, h = self._canvas_size()
        self._size = (w, h)
        rx, ry = w * scale / 2, h * scale / 2
        node_map = {}
        for i, v in enumerate(vertices):
            angle = 2 * math.pi * i / max(n, 1)
            x = w / 2 + (rx * math.cos(angle) if n > 1 else 0.0)
            y = h / 2 + (ry * math.sin(angle) if n > 1 else 0.0)
            node_map[v] = Node(self, x, y, labels[v] if labels else v)
        for u, v in G.edges():
            if u != v:
                Edge(self, node_map[u], node_map[v])
        return node_map

    def _on_resize(self, event):
        """Scale the drawing with the canvas so that the layout is kept."""
        if self._size and event.width > 1 and event.height > 1:
            sx = event.width / self._size[0]
            sy = event.height / self._size[1]
            for node in self.nodes:
                node.move_to(node.x * sx, node.y * sy)
        if event.width > 1 and event.height > 1:
            self._size = (event.width, event.height)


def _parse_graph6(text: str) -> nx.Graph:
    """Read a graph6 string typed by the user (header and whitespace allowed)."""
    text = text.strip()
    if text.startswith(">>graph6<<"):
        text = text[len(">>graph6<<"):]
    return nx.from_graph6_bytes(text.encode("ascii"))


def _vertex_labels(G: nx.Graph) -> dict:
    return {v: str(v) for v in G.nodes}


def _format_distribution(dist) -> str:
    if not dist:
        return "—"
    return "; ".join(f"{rank}: {count}" for rank, count in dist.items())


# ─────────────────────────────────────────────────────────────────────────────
#  Main editor
# ─────────────────────────────────────────────────────────────────────────────

class GraphEditor(ttk.Frame):
    """Graph editor with LC-orbit controls and a live statistics panel."""

    HELP_TEXT = (
        "Node mode: click the canvas to add a vertex, drag a vertex to move it,\n"
        "click a vertex or edge to select it, double-click a vertex to rename it.\n"
        "Edge mode: click two vertices to join them.\n"
        "Delete key or Delete button: remove the selected vertex or edge.\n\n"
        "Perform LC: complement the edges among the neighbours of the selected vertex.\n"
        "Compute Orbit: find every graph reachable by LC, up to isomorphism.\n"
        "Next / Prev in Orbit: step through the orbit.\n"
        "Orbit overview: draw the orbit as a graph; click a member to load it.\n"
        "Orbit stats: table of invariants for every member of the orbit."
    )

    def __init__(self, master):
        super().__init__(master)
        self.master = master
        master.title("Graph State Orbit Explorer")
        master.configure(bg=BG)
        self.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        self.edge_mode = False
        self.pending_node = None
        self.selected = None
        self._drag = None
        self.current_g6 = None
        self.orbit: list[str] = []
        self.orbit_neighbors: dict = {}
        self._orbit_index: dict[str, int] = {}
        self.overview = None
        self._orbit_running = False
        self._stats_running = False

        self._build_styles()
        self._build_toolbar()

        self.pane = ttk.Panedwindow(self, orient=tk.HORIZONTAL)
        self.pane.pack(fill=tk.BOTH, expand=True)
        self.view = GraphView(self.pane)
        self.right_frame = ttk.Frame(self.pane, width=300)
        self.pane.add(self.view, weight=4)
        self.pane.add(self.right_frame, weight=1)
        self._build_stats_panel(self.right_frame)

        canvas = self.view.canvas
        canvas.bind("<Button-1>", self._on_press)
        canvas.bind("<B1-Motion>", self._on_drag)
        canvas.bind("<ButtonRelease-1>", self._on_release)
        canvas.bind("<Double-Button-1>", self._on_double_click)
        master.bind("<Delete>", self._on_delete_key)

        self._graph_changed()

    # ----- UI construction -----
    def _build_styles(self):
        style = ttk.Style()
        try:
            style.theme_use("clam")
        except tk.TclError:
            pass
        style.configure("TFrame", background=BG)
        style.configure("TLabel", background=BG, foreground="white", font=("Helvetica", 11))
        style.configure("Small.TLabel", background=BG, foreground="#ddd", font=("Helvetica", 10))
        style.configure("Header.TLabel", background=BG, foreground="white",
                        font=("Helvetica", 12, "bold"))
        style.configure("TButton", padding=6)
        style.configure("Active.TButton", padding=6, background="#5050a0", foreground="white")
        style.map("Active.TButton", background=[("active", "#6868c0")])
        style.configure("Inactive.TButton", padding=6, background="#3a3a3a", foreground="white")
        style.map("Inactive.TButton", background=[("active", "#4a4a4a")])
        style.configure("TProgressbar", troughcolor="#1f1f1f", background="#68d391")

    def _build_toolbar(self):
        # Row 1: drawing modes and the LC orbit.
        bar = ttk.Frame(self)
        bar.pack(fill=tk.X, pady=(0, 4))
        self.toolbar = bar
        self.node_button = ttk.Button(bar, text="Node mode", style="Active.TButton",
                                      command=self.activate_node_mode)
        self.node_button.pack(side=tk.LEFT)
        self.edge_button = ttk.Button(bar, text="Edge mode", style="Inactive.TButton",
                                      command=self.activate_edge_mode)
        self.edge_button.pack(side=tk.LEFT, padx=(4, 0))
        ttk.Button(bar, text="Perform LC", command=self.perform_lc).pack(side=tk.LEFT, padx=(16, 0))
        self.orbit_button = ttk.Button(bar, text="Compute Orbit", command=self.compute_orbit)
        self.orbit_button.pack(side=tk.LEFT, padx=(4, 0))
        ttk.Button(bar, text="Prev in Orbit", command=lambda: self.step_orbit(-1)).pack(
            side=tk.LEFT, padx=(4, 0))
        ttk.Button(bar, text="Next in Orbit", command=lambda: self.step_orbit(1)).pack(side=tk.LEFT)
        ttk.Button(bar, text="Orbit overview", command=self.show_orbit_overview).pack(
            side=tk.LEFT, padx=(4, 0))
        self.stats_button = ttk.Button(bar, text="Orbit stats", command=self.show_orbit_stats)
        self.stats_button.pack(side=tk.LEFT, padx=(4, 0))

        # Row 2: editing, entanglement, colouring and help.
        bar2 = ttk.Frame(self)
        bar2.pack(fill=tk.X, pady=(0, 4))
        ttk.Button(bar2, text="Delete", command=self.delete_selected).pack(side=tk.LEFT)
        ttk.Button(bar2, text="Load graph6", command=self.load_graph6).pack(side=tk.LEFT, padx=(4, 0))
        ttk.Button(bar2, text="Compute Schmidt", command=self.compute_schmidt).pack(
            side=tk.LEFT, padx=(16, 0))
        ttk.Label(bar2, text="k-colour?").pack(side=tk.LEFT, padx=(16, 4))
        self.kcolor_var = tk.IntVar(value=3)
        ttk.Spinbox(bar2, from_=2, to=10, width=4, textvariable=self.kcolor_var).pack(side=tk.LEFT)
        ttk.Button(bar2, text="Check", command=self.check_k_colorable).pack(side=tk.LEFT, padx=(4, 0))
        ttk.Button(bar2, text="Greedy colour", command=self.apply_greedy_coloring).pack(
            side=tk.LEFT, padx=(4, 0))
        ttk.Button(bar2, text="Reset colours", command=self.reset_colors).pack(side=tk.LEFT, padx=(4, 0))
        ttk.Button(bar2, text="Help",
                   command=lambda: messagebox.showinfo("Help", self.HELP_TEXT)).pack(side=tk.RIGHT)

        # Status line.
        line = ttk.Frame(self)
        line.pack(fill=tk.X, pady=(0, 4))
        self.status_var = tk.StringVar(value="Mode: Node")
        ttk.Label(line, textvariable=self.status_var).pack(side=tk.LEFT, padx=(4, 16))
        self.graph6_var = tk.StringVar(value="Graph6: —")
        ttk.Label(line, textvariable=self.graph6_var, style="Small.TLabel").pack(side=tk.LEFT, padx=(0, 16))
        self.schmidt_var = tk.StringVar(value="Schmidt: —")
        ttk.Label(line, textvariable=self.schmidt_var, style="Small.TLabel").pack(side=tk.LEFT)

    def _attach_tooltip(self, widget, key: str):
        if not hasattr(self, "_tooltips"):
            self._tooltips = []
        self._tooltips.append(Tooltip(widget, STATS_HELP.get(key, key)))

    def _build_stats_panel(self, parent):
        """Right-hand panel: graph-state quantities first, then graph invariants. Scrolls."""
        canvas = tk.Canvas(parent, bg=BG, highlightthickness=0, width=300)
        scrollbar = ttk.Scrollbar(parent, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=scrollbar.set)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        panel = ttk.Frame(canvas)
        window = canvas.create_window((0, 0), window=panel, anchor="nw")
        panel.bind("<Configure>", lambda _e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.bind("<Configure>", lambda e: canvas.itemconfigure(window, width=e.width))
        self._bind_wheel(parent, canvas)

        ttk.Label(panel, text="Graph state", style="Header.TLabel").pack(anchor="w", padx=8, pady=(8, 4))
        self.orbit_size_var = tk.StringVar(value="Orbit size: —")
        orbit_label = ttk.Label(panel, textvariable=self.orbit_size_var, style="Small.TLabel",
                                wraplength=280)
        orbit_label.pack(anchor="w", padx=8, pady=(0, 4))
        self._attach_tooltip(orbit_label, "_orbit_size")

        self.schmidt_tv = ttk.Treeview(panel, columns=("rank", "count"), show="headings", height=5)
        self.schmidt_tv.heading("rank", text="Schmidt rank")
        self.schmidt_tv.heading("count", text="Bipartitions")
        self.schmidt_tv.column("rank", anchor="center", width=110)
        self.schmidt_tv.column("count", anchor="center", width=110)
        self.schmidt_tv.pack(fill=tk.X, padx=8)
        self._attach_tooltip(self.schmidt_tv, "_schmidt_tv")

        ttk.Label(panel, text="Graph stats", style="Header.TLabel").pack(anchor="w", padx=8, pady=(12, 4))
        grid = ttk.Frame(panel)
        grid.pack(fill=tk.BOTH, expand=True, padx=8)
        rows = [
            ("schmidt_rank_max", "Schmidt rank (max)"),
            ("nodes", "Nodes"), ("edges", "Edges"), ("density", "Density"),
            ("clustering", "Clustering"), ("assortativity", "Assortativity"),
            ("connected", "Connected"), ("components", "Components"),
            ("diameter_giant", "Diameter (giant)"), ("avg_path_len_giant", "Avg. path len (giant)"),
            ("bipartite", "Bipartite"), ("bipartition_sizes", "Bipartition sizes"),
            ("planar", "Planar"), ("clique_number", "Clique # (ω)"),
            ("chromatic_lower_bound", "χ lower bound"), ("chromatic_upper_bound", "χ upper bound"),
            ("deg_min", "deg min"), ("deg_avg", "deg avg"), ("deg_max", "deg max"),
            ("triangles", "Triangles"),
        ]
        self._stat_labels = {}
        self._density_bar = ttk.Progressbar(grid, orient="horizontal", mode="determinate",
                                            length=120, maximum=1.0)
        self._cluster_bar = ttk.Progressbar(grid, orient="horizontal", mode="determinate",
                                            length=120, maximum=1.0)
        for r, (key, text) in enumerate(rows):
            name = ttk.Label(grid, text=text)
            name.grid(row=r, column=0, sticky="w", pady=1)
            if key == "density":
                widget = self._density_bar
                widget.grid(row=r, column=1, sticky="ew", pady=1)
            elif key == "clustering":
                widget = self._cluster_bar
                widget.grid(row=r, column=1, sticky="ew", pady=1)
            else:
                widget = ttk.Label(grid, text="—", style="Small.TLabel")
                widget.grid(row=r, column=1, sticky="e", pady=1)
            self._stat_labels[key] = widget
            self._attach_tooltip(name, key)
            self._attach_tooltip(widget, key)
        grid.columnconfigure(1, weight=1)

        buttons = ttk.Frame(panel)
        buttons.pack(fill=tk.X, padx=8, pady=8)
        ttk.Button(buttons, text="Refresh stats", command=self._graph_changed).pack(side=tk.LEFT)
        ttk.Button(buttons, text="Copy stats CSV", command=self.copy_stats_csv).pack(
            side=tk.LEFT, padx=6)

    @staticmethod
    def _bind_wheel(area, canvas):
        """Scroll ``canvas`` with the mouse wheel while the pointer is over ``area``."""
        def on_wheel(event):
            if event.num == 4 or event.delta > 0:
                canvas.yview_scroll(-1, "units")
            elif event.num == 5 or event.delta < 0:
                canvas.yview_scroll(1, "units")

        def enter(_e):
            for sequence in ("<MouseWheel>", "<Button-4>", "<Button-5>"):
                area.bind_all(sequence, on_wheel)

        def leave(_e):
            inside = area.winfo_containing(*area.winfo_pointerxy())
            if inside is not None and str(inside).startswith(str(area)):
                return  # the pointer moved onto a child widget of the panel
            for sequence in ("<MouseWheel>", "<Button-4>", "<Button-5>"):
                area.unbind_all(sequence)

        area.bind("<Enter>", enter, add="+")
        area.bind("<Leave>", leave, add="+")

    # ----- modes and selection -----
    def activate_node_mode(self):
        self.edge_mode = False
        self._reset_pending()
        self.status_var.set("Mode: Node")
        self.node_button.config(style="Active.TButton")
        self.edge_button.config(style="Inactive.TButton")

    def activate_edge_mode(self):
        self.edge_mode = True
        self._select(None)
        self.status_var.set("Mode: Edge – pick two vertices")
        self.edge_button.config(style="Active.TButton")
        self.node_button.config(style="Inactive.TButton")

    def _select(self, item):
        if self.selected is not None and self.selected is not item:
            self.selected.set_selected(False)
        self.selected = item
        if item is not None:
            item.set_selected(True)

    def _reset_pending(self):
        if self.pending_node is not None:
            self.pending_node.set_selected(False)
        self.pending_node = None

    # ----- mouse and keyboard -----
    def _on_press(self, event):
        item = self.view.item_at(event.x, event.y)
        if self.edge_mode:
            if isinstance(item, Node):
                self._edge_click(item)
            else:
                self._reset_pending()
                self.status_var.set("Mode: Edge – pick two vertices")
            return
        if isinstance(item, Node):
            self._select(item)
            self._drag = (item, event.x, event.y)
        elif isinstance(item, Edge):
            self._select(item)
        else:
            self._select(None)
            self._add_node(event.x, event.y)

    def _on_drag(self, event):
        if self._drag is None:
            return
        node, x0, y0 = self._drag
        node.move_to(node.x + event.x - x0, node.y + event.y - y0)
        self._drag = (node, event.x, event.y)

    def _on_release(self, _event):
        self._drag = None

    def _on_double_click(self, event):
        if self.edge_mode:
            # A quick second click, e.g. to start a new edge at the vertex that
            # just ended one, is an ordinary click in edge mode.
            self._on_press(event)
            return
        node = self.view.item_at(event.x, event.y)
        if not isinstance(node, Node):
            return
        new = simpledialog.askstring("Rename vertex", "New label:", initialvalue=node.label,
                                     parent=self.master)
        if new is None or not new.strip():
            return
        new = new.strip()
        if any(other is not node and other.label == new for other in self.view.nodes):
            messagebox.showwarning("Rename vertex", f"Another vertex is already labelled {new!r}.")
            return
        node.set_label(new)

    def _on_delete_key(self, event):
        if isinstance(event.widget, (tk.Entry, ttk.Entry)):
            return
        self.delete_selected()

    def _edge_click(self, node: Node):
        if self.pending_node is None:
            self.pending_node = node
            node.set_selected(True)
            self.status_var.set("Mode: Edge – pick second vertex")
            return
        first = self.pending_node
        self._reset_pending()
        self.status_var.set("Mode: Edge – pick two vertices")
        try:
            Edge(self.view, first, node)
        except ValueError as exc:
            messagebox.showwarning("Invalid edge", str(exc))
            return
        self._graph_changed()

    # ----- editing -----
    def _next_label(self) -> str:
        used = {node.label for node in self.view.nodes}
        numbers = [int(label) for label in used if label.isdigit()]
        candidate = max(numbers) + 1 if numbers else 0
        while str(candidate) in used:
            candidate += 1
        return str(candidate)

    def _add_node(self, x: float, y: float):
        Node(self.view, x, y, self._next_label())
        self._graph_changed()

    def delete_selected(self):
        item = self.selected
        if item is None:
            return
        self._select(None)
        if item is self.pending_node:
            self._reset_pending()
        item.delete()
        self._graph_changed()

    def show_graph(self, G: nx.Graph):
        """Draw ``G`` on a circle, replacing the current drawing."""
        self._select(None)
        self._reset_pending()
        self.view.draw(G, labels=_vertex_labels(G))
        self._graph_changed()

    def perform_lc(self):
        """Local complementation at the selected vertex, edited in place."""
        node = self.selected if isinstance(self.selected, Node) else None
        if node is None:
            messagebox.showinfo("Local complementation",
                                "Select a vertex first (click it in node mode).")
            return
        neighbours = [edge.other(node) for edge in node.edges]
        for a, b in combinations(neighbours, 2):
            edge = self.view.edge_between(a, b)
            if edge is not None:
                edge.delete()
            else:
                Edge(self.view, a, b)
        self._graph_changed()

    # ----- state refresh -----
    def _graph_changed(self):
        """Recompute everything that depends on the drawn graph."""
        G = self.view.to_networkx()
        self.current_g6 = canonical_relabel(G).decode("ascii") if G.number_of_nodes() else None
        self.graph6_var.set(f"Graph6: {self.current_g6 or '—'}")
        self._update_stats_panel(G)
        self._update_orbit_label()
        if self.overview is not None:
            self.overview.follow(self.current_g6)

    def _update_orbit_label(self):
        if not self.orbit:
            self.orbit_size_var.set("Orbit size: —")
        elif self.current_g6 in self._orbit_index:
            self.orbit_size_var.set(f"Orbit size: {len(self.orbit)}")
        else:
            self.orbit_size_var.set(f"Orbit size: {len(self.orbit)} (current graph not in this orbit)")

    def _update_stats_panel(self, G: nx.Graph):
        s = graph_stats(G, schmidt_max_nodes=LIVE_SCHMIDT_MAX_NODES)
        self._last_stats = s
        for key, widget in self._stat_labels.items():
            if key == "density":
                widget["value"] = s["density"]
            elif key == "clustering":
                widget["value"] = s["clustering"]
            else:
                widget.config(text=self._format_stat(key, s.get(key)))
        if s["schmidt_rank_max"] is not None:
            self._show_schmidt(s["schmidt_rank_max"], s["schmidt_distribution"])
        else:
            self.schmidt_tv.delete(*self.schmidt_tv.get_children())
            self.schmidt_var.set("Schmidt: —" if s["nodes"] == 0
                                 else "Schmidt: press Compute Schmidt")

    @staticmethod
    def _format_stat(key, value):
        if value is None:
            return "—"
        if key == "bipartition_sizes":
            return f"{value[0]} | {value[1]}"
        if isinstance(value, float):
            return f"{value:.4f}" if abs(value) < 10 else f"{value:.2f}"
        return str(value)

    def _show_schmidt(self, max_rank: int, dist: dict):
        self.schmidt_var.set(f"Schmidt rank (max): {max_rank}")
        self.schmidt_tv.delete(*self.schmidt_tv.get_children())
        for rank, count in dist.items():
            self.schmidt_tv.insert("", "end", values=(rank, count))

    # ----- orbit -----
    def compute_orbit(self):
        """Explore the LC orbit of the current graph in a worker thread."""
        if self._orbit_running:
            return
        G = self.view.to_networkx()
        if G.number_of_nodes() == 0:
            messagebox.showinfo("Orbit", "Draw or load a graph first.")
            return
        seed = nx.to_graph6_bytes(G, header=False)
        self._orbit_running = True
        self.orbit_button.state(["disabled"])

        def worker():
            try:
                orbit, neighbors = graphOperations.exploreOrbit([seed], self.master)
            except Exception as exc:
                message = str(exc)
                self._post(lambda: self._orbit_failed(message))
                return
            self._post(lambda: self._orbit_done(orbit, neighbors))

        threading.Thread(target=worker, daemon=True).start()

    def _post(self, callback):
        """Schedule ``callback`` on the Tk main thread (safe from worker threads)."""
        try:
            self.master.after(0, callback)
        except (RuntimeError, tk.TclError):
            pass  # the window has been closed

    def _orbit_failed(self, message: str):
        self._orbit_running = False
        self.orbit_button.state(["!disabled"])
        messagebox.showerror("Orbit error", message)

    def _orbit_done(self, orbit: list, neighbors: dict):
        self._orbit_running = False
        self.orbit_button.state(["!disabled"])
        self.orbit = orbit
        self.orbit_neighbors = neighbors
        self._orbit_index = {g6: i for i, g6 in enumerate(orbit)}
        if self.overview is not None:
            self.overview.close()
        self._update_orbit_label()

    def step_orbit(self, step: int):
        """Show the next (step=1) or previous (step=-1) member of the orbit."""
        if not self.orbit:
            messagebox.showinfo("Orbit", "Compute an orbit first.")
            return
        if self.current_g6 in self._orbit_index:
            i = (self._orbit_index[self.current_g6] + step) % len(self.orbit)
        else:
            i = 0 if step > 0 else len(self.orbit) - 1
        self.show_graph(nx.from_graph6_bytes(self.orbit[i].encode("ascii")))

    def show_orbit_overview(self):
        if not self.orbit:
            messagebox.showinfo("No orbit", "Compute an orbit first (Compute Orbit button).")
            return
        if self.overview is not None:
            self.overview.lift()
            return
        self.overview = OrbitOverview(self)

    # ----- orbit statistics -----
    def show_orbit_stats(self):
        """Compute the statistics table of the orbit in the background, then show it."""
        if not self.orbit:
            messagebox.showinfo("No orbit", "Compute an orbit first.")
            return
        if self._stats_running:
            return
        graphs = list(self.orbit)
        self._stats_running = True
        self.stats_button.state(["disabled"])
        self.status_var.set(f"Computing statistics for {len(graphs)} graphs…")

        def worker():
            try:
                if len(graphs) >= PARALLEL_MIN_ROWS:
                    context = multiprocessing.get_context("spawn")
                    with ProcessPoolExecutor(mp_context=context) as pool:
                        rows = list(pool.map(orbit_table_row, graphs, chunksize=32))
                else:
                    rows = [orbit_table_row(g6) for g6 in graphs]
            except Exception as exc:
                message = str(exc)
                self._post(lambda: self._orbit_stats_finished(None, message))
                return
            self._post(lambda: self._orbit_stats_finished(rows, None))

        threading.Thread(target=worker, daemon=True).start()

    def _orbit_stats_finished(self, rows, error):
        self._stats_running = False
        self.stats_button.state(["!disabled"])
        self.status_var.set("Mode: Edge – pick two vertices" if self.edge_mode else "Mode: Node")
        if error is not None:
            messagebox.showerror("Orbit stats", error)
            return
        self._open_orbit_stats_window(rows)

    def _open_orbit_stats_window(self, rows):
        top = tk.Toplevel(self)
        top.title(f"Orbit stats ({len(rows)} graphs)")
        tree = ttk.Treeview(top, columns=ORBIT_TABLE_COLUMNS, show="headings", height=18)
        tree.pack(fill=tk.BOTH, expand=True)
        for column in ORBIT_TABLE_COLUMNS:
            tree.heading(column, text=column)
            tree.column(column, anchor="center", width=280 if column == "graph6" else 110)
        for row in rows:
            tree.insert("", "end", values=["—" if v is None else v for v in row])
        add_treeview_header_tooltips(tree, ORBIT_STATS_COL_HELP)

        def export_csv():
            path = filedialog.asksaveasfilename(
                parent=top, title="Export orbit statistics", defaultextension=".csv",
                filetypes=[("CSV files", "*.csv")], initialfile="orbit_stats.csv")
            if not path:
                return
            with open(path, "w", newline="", encoding="utf-8") as f:
                writer = csv.writer(f)
                writer.writerow(ORBIT_TABLE_COLUMNS)
                writer.writerows(rows)
            messagebox.showinfo("Export", f"Saved {path}", parent=top)

        ttk.Button(top, text="Export CSV", command=export_csv).pack(pady=6)

    # ----- dialogs and utilities -----
    def load_graph6(self):
        """Ask for a graph6 string and draw the graph."""
        text = simpledialog.askstring("Load graph6", "graph6 string:", parent=self.master)
        if not text:
            return
        try:
            G = _parse_graph6(text)
        except Exception as exc:
            messagebox.showerror("Load error", f"Not a valid graph6 string: {exc}")
            return
        self.show_graph(G)

    def compute_schmidt(self):
        """Enumerate all bipartitions of the current graph state and report Schmidt ranks."""
        G = self.view.to_networkx()
        n = G.number_of_nodes()
        if n == 0:
            return
        if n > BUTTON_SCHMIDT_MAX_NODES:
            messagebox.showwarning(
                "Schmidt rank",
                f"Exhaustive enumeration covers 2^(n−1) − 1 bipartitions and is limited "
                f"to n ≤ {BUTTON_SCHMIDT_MAX_NODES} vertices (this graph has {n}).")
            return
        self.master.config(cursor="watch")
        self.master.update_idletasks()
        try:
            max_rank, dist = schmidt_rank_profile(G)
        finally:
            self.master.config(cursor="")
        self._show_schmidt(max_rank, dist)
        ebits = int(math.log2(max_rank))
        messagebox.showinfo(
            "Schmidt rank",
            f"Maximum Schmidt rank: {max_rank} ({ebits} ebits)\n"
            f"Bipartitions: {sum(dist.values())}\n"
            f"Schmidt rank: count  →  {_format_distribution(dist)}")

    def _apply_coloring(self, coloring: dict):
        for i, node in enumerate(self.view.nodes):
            if i in coloring:
                node.set_fill(PALETTE[coloring[i] % len(PALETTE)])

    def apply_greedy_coloring(self):
        self._apply_coloring(greedy_coloring(self.view.to_networkx()))

    def reset_colors(self):
        for node in self.view.nodes:
            node.set_fill(Node.FILL)

    def check_k_colorable(self):
        """Exact k-colourability test; on success the colouring is drawn."""
        try:
            k = int(self.kcolor_var.get())
        except (tk.TclError, ValueError):
            messagebox.showwarning("k-colourable", "Enter an integer k.")
            return
        status, coloring = is_k_colorable(self.view.to_networkx(), k)
        if status is True:
            self._apply_coloring(coloring)
            messagebox.showinfo("k-colourable", f"Yes, the graph is {k}-colourable. "
                                                "The colouring is shown.")
        elif status is False:
            messagebox.showinfo("k-colourable", f"No, the graph is not {k}-colourable.")
        else:
            messagebox.showwarning("k-colourable", "Undecided: the exact search stopped "
                                                   "at its step limit.")

    def copy_stats_csv(self):
        """Copy the current statistics to the clipboard as CSV."""
        s = self._last_stats
        buffer = io.StringIO()
        writer = csv.writer(buffer, lineterminator="\n")
        writer.writerow(("metric", "value"))
        for key, value in s.items():
            if key == "schmidt_distribution":
                value = _format_distribution(value)
            elif isinstance(value, tuple):
                value = f"{value[0]}|{value[1]}"
            writer.writerow((key, "" if value is None else value))
        text = buffer.getvalue()
        self.clipboard_clear()
        self.clipboard_append(text)
        messagebox.showinfo("Copied", "Stats copied to the clipboard as CSV.")


# ─────────────────────────────────────────────────────────────────────────────
#  Orbit overview
# ─────────────────────────────────────────────────────────────────────────────

class OrbitOverview(tk.Toplevel):
    """The LC orbit drawn as a graph: vertices are orbit members, edges single LC moves.

    LC moves that map a graph to an isomorphic copy of itself are not drawn.
    Click a member to load it in the editor and show its neighbourhood up to
    the chosen depth; click the empty canvas to show the whole orbit.
    """

    def __init__(self, editor: GraphEditor):
        super().__init__(editor.master)
        self.editor = editor
        self.title(f"Orbit overview ({len(editor.orbit)} graphs)")
        self.geometry("800x600")
        self.configure(bg=BG)

        self.meta = nx.Graph()
        self.meta.add_nodes_from(editor.orbit)
        for source, moves in editor.orbit_neighbors.items():
            for target, _vertex in moves:
                if target != source:
                    self.meta.add_edge(source, target)

        controls = ttk.Frame(self)
        controls.pack(fill=tk.X)
        ttk.Label(controls, text="Depth").pack(side=tk.LEFT, padx=8, pady=6)
        self.depth_var = tk.IntVar(value=2)
        ttk.Spinbox(controls, from_=1, to=10, width=4, textvariable=self.depth_var).pack(side=tk.LEFT)
        self.info_var = tk.StringVar()
        ttk.Label(controls, textvariable=self.info_var, style="Small.TLabel").pack(
            side=tk.LEFT, padx=12)

        self.view = GraphView(self)
        self.view.pack(fill=tk.BOTH, expand=True)
        self.view.canvas.bind("<Button-1>", self._on_press)
        self.view.canvas.bind("<B1-Motion>", self._on_drag)
        self.view.canvas.bind("<ButtonRelease-1>", self._on_release)
        self._press = None
        self.protocol("WM_DELETE_WINDOW", self.close)

        self.update_idletasks()
        start = None
        if self.meta.number_of_nodes() > OVERVIEW_FULL_MAX_NODES:
            start = editor.current_g6 if editor.current_g6 in self.meta else editor.orbit[0]
        self.show(start)

    def _depth(self) -> int:
        try:
            return max(1, int(self.depth_var.get()))
        except (tk.TclError, ValueError):
            return 1

    def show(self, center=None):
        """Draw the whole orbit (center=None) or the neighbourhood of ``center``."""
        if center is None:
            if self.meta.number_of_nodes() > OVERVIEW_FULL_MAX_NODES:
                self.info_var.set(f"Orbit has more than {OVERVIEW_FULL_MAX_NODES} members; "
                                  "click a member to see its neighbourhood.")
                return
            nodes = list(self.meta.nodes)
            self.info_var.set("Whole orbit – click a member for its neighbourhood")
        else:
            depths = nx.single_source_shortest_path_length(self.meta, center, cutoff=self._depth())
            nodes = [v for v in self.meta.nodes if v in depths]
            self.info_var.set(f"Neighbourhood of depth {self._depth()} – click the canvas "
                              "for the whole orbit")
        self.view.draw(self.meta.subgraph(nodes), labels={v: v for v in nodes})
        self.highlight(self.editor.current_g6)

    def highlight(self, g6):
        for node in self.view.nodes:
            node.set_marked(node.label == g6)

    def follow(self, g6):
        """Highlight the editor's current graph, redrawing around it if it is hidden."""
        if g6 in self.meta and all(node.label != g6 for node in self.view.nodes):
            self.show(g6)
        else:
            self.highlight(g6)

    def _on_press(self, event):
        item = self.view.item_at(event.x, event.y)
        node = item if isinstance(item, Node) else None
        self._press = [node, event.x, event.y, False]

    def _on_drag(self, event):
        if not self._press or self._press[0] is None:
            return
        node, x0, y0, _moved = self._press
        node.move_to(node.x + event.x - x0, node.y + event.y - y0)
        self._press = [node, event.x, event.y, True]

    def _on_release(self, _event):
        if not self._press:
            return
        node, _x, _y, moved = self._press
        self._press = None
        if moved:
            return
        if node is None:
            self.show(None)
            return
        g6 = node.label
        self.show(g6)
        self.editor.show_graph(nx.from_graph6_bytes(g6.encode("ascii")))

    def close(self):
        self.editor.overview = None
        self.destroy()


# ─────────────────────────────────────────────────────────────────────────────
#  Entry point
# ─────────────────────────────────────────────────────────────────────────────

def main():
    root = tk.Tk()
    root.geometry("1100x700")
    root.minsize(800, 500)
    editor = GraphEditor(root)
    editor.show_graph(nx.cycle_graph(5))  # five-qubit ring state as a first example
    root.mainloop()


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
