# -*- coding: utf-8 -*-
"""
TF-Explorer Interactive SVG Genome Browser
Generates ultra-premium, zero-dependency interactive vector graphics with hover tooltips for promoter tracks.
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Optional

def generate_interactive_genome_svg(
    df_encode: pd.DataFrame,
    df_motifs: pd.DataFrame,
    df_cons: Optional[pd.DataFrame],
    df_synergy: Optional[pd.DataFrame],
    promoter_up: int,
    promoter_down: int,
    gene_name: str,
    top_n: Optional[int] = None,
    threshold: Optional[float] = None
) -> str:
    """
    Constructs a responsive HTML string containing an interactive SVG and inline JavaScript
    for premium mouse hover tooltips and micro-animations.
    """
    
    # 1. Dimensions and Grid Setup
    svg_w = 1000
    svg_h = 390
    margin_l = 70
    margin_r = 70
    track_w = svg_w - margin_l - margin_r  # 860px
    total_bp = promoter_up + promoter_down
    
    def get_x(bp_offset: float) -> float:
        # bp_offset is distance to TSS (-promoter_up to +promoter_down)
        rel = bp_offset + promoter_up
        return margin_l + (rel * track_w / total_bp)
        
    tss_x = get_x(0)
    
    # 2. Setup CSS Styles and Theme (Dark Glassmorphic Theme)
    svg_styles = """
    <style>
        @import url('https://fonts.googleapis.com/css2?family=Outfit:wght@300;400;500;600;700&family=Inter:wght@300;400;500;600&display=swap');
        
        .svg-card {
            font-family: 'Outfit', 'Inter', -apple-system, sans-serif;
            background: #0f141c;
            border-radius: 14px;
            padding: 20px;
            border: 1px solid rgba(255, 255, 255, 0.08);
            box-shadow: 0 12px 40px rgba(0, 0, 0, 0.4);
            position: relative;
            user-select: none;
        }
        
        .svg-title {
            color: #f1f5f9;
            font-size: 16px;
            font-weight: 600;
            margin-bottom: 5px;
            letter-spacing: 0.5px;
            display: flex;
            align-items: center;
            gap: 8px;
        }
        
        .svg-subtitle {
            color: #94a3b8;
            font-size: 12px;
            margin-bottom: 20px;
        }
        
        .track-bg {
            fill: #151c27;
            rx: 6px;
        }
        
        .tss-line {
            stroke: #ef4444;
            stroke-dasharray: 4 4;
            stroke-width: 1.5;
        }
        
        .axis-label {
            fill: #64748b;
            font-size: 10px;
            font-weight: 500;
            font-family: 'Inter', sans-serif;
        }
        
        .track-header {
            fill: #94a3b8;
            font-size: 11px;
            font-weight: 600;
            letter-spacing: 0.5px;
        }
        
        /* Motif Hover Micro-animations */
        .motif-pill {
            cursor: pointer;
            rx: 4px;
            ry: 4px;
            transition: all 0.2s ease;
            fill: #3b82f6;
            fill-opacity: 0.6;
            stroke: #60a5fa;
            stroke-width: 1;
        }
        .motif-pill:hover {
            fill-opacity: 0.9;
            stroke: #93c5fd;
            stroke-width: 1.5;
            filter: drop-shadow(0 0 4px rgba(96, 165, 250, 0.6));
        }
        
        /* Peak Triangles */
        .peak-triangle {
            cursor: pointer;
            fill: #f43f5e;
            fill-opacity: 0.7;
            stroke: #fda4af;
            stroke-width: 1;
            transition: all 0.2s ease;
        }
        .peak-triangle:hover {
            fill-opacity: 1;
            stroke: #fecdd3;
            stroke-width: 1.5;
            transform: scale(1.2);
            transform-origin: center;
            filter: drop-shadow(0 0 6px rgba(244, 63, 94, 0.8));
        }
        .peak-triangle-high {
            fill: #ef4444;
            fill-opacity: 1;
            stroke: #ffffff;
            stroke-width: 1.2;
        }
        
        /* Synergy Region Glow */
        .synergy-band {
            fill: #f59e0b;
            fill-opacity: 0.12;
            stroke: #f59e0b;
            stroke-dasharray: 2 2;
            stroke-width: 0.8;
            stroke-opacity: 0.4;
            transition: fill-opacity 0.3s ease;
        }
        .synergy-band:hover {
            fill-opacity: 0.22;
        }
        
        /* Tooltip style */
        .custom-tooltip {
            position: absolute;
            background: rgba(15, 20, 28, 0.96);
            color: #f8fafc;
            padding: 10px 14px;
            border-radius: 8px;
            border: 1px solid rgba(255, 255, 255, 0.15);
            font-size: 11px;
            font-family: 'Inter', sans-serif;
            box-shadow: 0 8px 30px rgba(0, 0, 0, 0.6);
            pointer-events: none;
            backdrop-filter: blur(8px);
            z-index: 9999;
            transition: opacity 0.1s ease, transform 0.1s ease;
            opacity: 0;
            display: none;
            line-height: 1.5;
        }
        
        .tooltip-title {
            font-weight: 700;
            color: #60a5fa;
            margin-bottom: 4px;
            font-size: 12px;
        }
        .tooltip-row {
            margin-bottom: 2px;
        }
        .tooltip-label {
            color: #94a3b8;
            font-weight: 500;
        }
        .tooltip-val {
            color: #f1f5f9;
            font-weight: 600;
        }
    </style>
    """
    
    # Start HTML wrapper
    html_out = []
    html_out.append(f'<div class="svg-card">')
    html_out.append(f'<div class="svg-title">🧬 {gene_name} Interactive Promoter Track & Synergy Browser</div>')
    html_out.append(f'<div class="svg-subtitle">Hover over motifs, conservation waves, or ChIP peaks to query specific values. Grid shows promoter span from -{promoter_up} bp to +{promoter_down} bp relative to TSS.</div>')
    
    # Main SVG
    html_out.append(f'<svg viewBox="0 0 {svg_w} {svg_h}" style="width:100%; height:auto; overflow:visible;">')
    
    # Inject styling
    html_out.append(svg_styles)
    
    # 3. Definitions (Gradients and drop shadows)
    html_out.append("""
    <defs>
        <!-- Conservation Gradient -->
        <linearGradient id="consGrad" x1="0" y1="0" x2="0" y2="1">
            <stop offset="0%" stop-color="#14b8a6" stop-opacity="0.5"/>
            <stop offset="100%" stop-color="#14b8a6" stop-opacity="0.02"/>
        </linearGradient>
        <!-- Synergy Gradient -->
        <linearGradient id="synGrad" x1="0" y1="0" x2="0" y2="1">
            <stop offset="0%" stop-color="#f59e0b" stop-opacity="0.25"/>
            <stop offset="100%" stop-color="#ef4444" stop-opacity="0.02"/>
        </linearGradient>
    </defs>
    """)
    
    # 4. Draw Layout Grid Tracks
    # Track coordinates mapping
    # Track 1: Conservation (Y=85 to 165)
    # Track 2: TF Motifs (Y=185 to 245)
    # Track 3: ChIP Peaks (Y=265 to 325)
    html_out.append(f'<rect x="{margin_l}" y="85" width="{track_w}" height="80" class="track-bg" />')
    html_out.append(f'<rect x="{margin_l}" y="185" width="{track_w}" height="60" class="track-bg" />')
    html_out.append(f'<rect x="{margin_l}" y="265" width="{track_w}" height="60" class="track-bg" />')
    
    # Track Header labels (Left side of track boxes)
    html_out.append(f'<text x="15" y="130" class="track-header" fill="#14b8a6" style="fill: #14b8a6;">CONS</text>')
    html_out.append(f'<text x="15" y="220" class="track-header" fill="#3b82f6" style="fill: #3b82f6;">MOTIFS</text>')
    html_out.append(f'<text x="15" y="300" class="track-header" fill="#f43f5e" style="fill: #f43f5e;">PEAKS</text>')
    
    # 5. Plot Synergy Hotspots (drawn in background across all tracks)
    if df_synergy is not None and not df_synergy.empty:
        for _, row in df_synergy.iterrows():
            start_rel = row['start'] - promoter_up
            end_rel = row['end'] - promoter_up
            
            x_start = get_x(start_rel)
            x_end = get_x(end_rel)
            span_w = max(4.0, x_end - x_start)
            
            tooltip_txt = (
                f'<div class="tooltip-title" style="color: #f59e0b;">🌟 Synergy Hotspot</div>'
                f'<div class="tooltip-row"><span class="tooltip-label">Region:</span> <span class="tooltip-val">{start_rel} to {end_rel} bp</span></div>'
                f'<div class="tooltip-row"><span class="tooltip-label">TFs:</span> <span class="tooltip-val">{row["tfs"]}</span></div>'
                f'<div class="tooltip-row"><span class="tooltip-label">Motifs:</span> <span class="tooltip-val">{row["motifs_count"]} sites</span></div>'
            )
            
            html_out.append(f'<rect x="{x_start}" y="85" width="{span_w}" height="240" fill="url(#synGrad)" class="synergy-band" '
                            f'onmousemove="showTooltip(event, \'{tooltip_txt}\')" onmouseout="hideTooltip()" />')
            
            # Synergy horizontal indicator line at bottom
            html_out.append(f'<line x1="{x_start}" y1="327" x2="{x_end}" y2="327" stroke="#f59e0b" stroke-width="2.5" />')
            html_out.append(f'<text x="{(x_start+x_end)/2}" y="337" class="axis-label" style="fill:#f59e0b; font-weight:700;" text-anchor="middle">Synergy</text>')
            
    # 6. Plot Conservation Wave Track
    # Base-by-base phastCons score line and shaded path
    has_cons = df_cons is not None and not df_cons.empty
    if has_cons:
        df_cons_sorted = df_cons.sort_values("distance_to_tss")
        x_coords = []
        y_coords = []
        
        # We need to construct points for the polygon (shaded under-curve) and polyline (curve border)
        for _, row in df_cons_sorted.iterrows():
            rel = row['distance_to_tss']
            val = row['phastCons']
            if pd.isna(val):
                val = 0.0
            
            cx = get_x(rel)
            # Track goes from Y=85 (phastCons=1.0) to Y=165 (phastCons=0.0)
            cy = 165 - (val * 75)
            
            x_coords.append(cx)
            y_coords.append(cy)
            
        if x_coords:
            # Filled area
            points_fill = [f"{margin_l},165"]
            for cx, cy in zip(x_coords, y_coords):
                points_fill.append(f"{cx:.1f},{cy:.1f}")
            points_fill.append(f"{margin_l + track_w},165")
            
            html_out.append(f'<polygon points="{" ".join(points_fill)}" fill="url(#consGrad)" />')
            
            # Top stroke line
            points_stroke = []
            for cx, cy in zip(x_coords, y_coords):
                points_stroke.append(f"{cx:.1f},{cy:.1f}")
            html_out.append(f'<polyline points="{" ".join(points_stroke)}" fill="none" stroke="#14b8a6" stroke-width="1.8" />')
            
            # A dynamic mouse-tracking transparent layer to scan conservation base values!
            # We can sample every 5th index or build invisible interactive vertical strips for tooltips
            step_count = min(150, len(df_cons_sorted))
            indices = np.linspace(0, len(df_cons_sorted) - 1, step_count, dtype=int)
            for idx in indices:
                row = df_cons_sorted.iloc[idx]
                rel = row['distance_to_tss']
                phast = row['phastCons']
                phyl = row['phyloP']
                
                cx = get_x(rel)
                cx_w = max(2.0, track_w / step_count)
                
                tooltip_txt = (
                    f'<div class="tooltip-title" style="color: #14b8a6;">🧬 Conservation Score</div>'
                    f'<div class="tooltip-row"><span class="tooltip-label">Distance:</span> <span class="tooltip-val">{rel:+,} bp</span></div>'
                    f'<div class="tooltip-row"><span class="tooltip-label">phastCons (100-way):</span> <span class="tooltip-val">{phast:.4f}</span></div>'
                    f'<div class="tooltip-row"><span class="tooltip-label">phyloP (100-way):</span> <span class="tooltip-val">{phyl:.4f}</span></div>'
                )
                html_out.append(f'<rect x="{cx - cx_w/2}" y="85" width="{cx_w}" height="80" fill="transparent" style="cursor:crosshair;" '
                                f'onmousemove="showTooltip(event, \'{tooltip_txt}\')" onmouseout="hideTooltip()" />')
                
    # 7. Draw Motifs Track (Pills/Capsules)
    if not df_motifs.empty:
        # Avoid overlapping pills vertically by stacking them if they overlap on x-axis
        lanes = {}  # Map index to lane y-offset
        
        # Sort motifs by length to handle them properly
        df_sorted = df_motifs.sort_values("start").copy()
        
        for idx, row in df_sorted.iterrows():
            start_rel = row['start'] - promoter_up
            end_rel = row['end'] - promoter_up
            
            x_start = get_x(start_rel)
            x_end = get_x(end_rel)
            w = max(5.0, x_end - x_start)
            
            # Simple lane packing algorithm to prevent overlaps
            # TF motif tracks Y range: 185 to 245
            lane = 0
            while True:
                # check if lane intersects with any pill in this lane
                collision = False
                if lane in lanes:
                    for prev_start, prev_end in lanes[lane]:
                        if not (x_end < prev_start or x_start > prev_end):
                            collision = True
                            break
                if not collision:
                    break
                lane += 1
                if lane >= 3:  # max 3 lanes to fit inside height 60
                    lane = 0
                    break
            
            if lane not in lanes:
                lanes[lane] = []
            lanes[lane].append((x_start, x_end))
            
            py = 190 + (lane * 18)
            ph = 14
            
            # Map different TFs to beautiful distinct colors
            tf_colors = {
                "YY1": ("#a78bfa", "#c4b5fd"),
                "CREB": ("#f472b6", "#fbcfe8"),
                "E2F1": ("#60a5fa", "#93c5fd"),
                "GAPDH": ("#34d399", "#6ee7b7"),
                "MYC": ("#fbbf24", "#fcd34d"),
                "SP1": ("#2dd4bf", "#99f6e4"),
            }
            
            color_pair = tf_colors.get(row['tf_name'], ("#3b82f6", "#60a5fa"))
            fill_c = color_pair[0]
            stroke_c = color_pair[1]
            
            tooltip_txt = (
                f'<div class="tooltip-title" style="color: {stroke_c};">🏷️ JASPAR TF Binding Motif</div>'
                f'<div class="tooltip-row"><span class="tooltip-label">TF Name:</span> <span class="tooltip-val">{row["tf_name"]}</span></div>'
                f'<div class="tooltip-row"><span class="tooltip-label">Motif ID:</span> <span class="tooltip-val">{row.get("jaspar_id", "N/A")}</span></div>'
                f'<div class="tooltip-row"><span class="tooltip-label">Sequence:</span> <span class="tooltip-val" style="font-family:monospace; color:#34d399; font-size:11px;">{row.get("sequence", "N/A")}</span></div>'
                f'<div class="tooltip-row"><span class="tooltip-label">PWM Score:</span> <span class="tooltip-val">{row.get("score", 0.0):.2f}</span></div>'
                f'<div class="tooltip-row"><span class="tooltip-label">Position:</span> <span class="tooltip-val">{start_rel:+} to {end_rel:+} bp</span></div>'
                f'<div class="tooltip-row"><span class="tooltip-label">Strand:</span> <span class="tooltip-val">{row.get("strand", "+")}</span></div>'
            )
            
            # Pill Rect
            html_out.append(f'<rect x="{x_start}" y="{py}" width="{w}" height="{ph}" class="motif-pill" '
                            f'style="fill: {fill_c}; stroke: {stroke_c};" '
                            f'onmousemove="showTooltip(event, \'{tooltip_txt}\')" onmouseout="hideTooltip()" />')
            
            # Centered Text in pill if width is enough
            if w > 25:
                html_out.append(f'<text x="{x_start + w/2}" y="{py + 10}" fill="#0f141c" font-size="9px" font-weight="700" '
                                f'text-anchor="middle" pointer-events="none">{row["tf_name"]}</text>')
                
    # 8. Draw ChIP Peaks Track
    if not df_encode.empty:
        window_peaks = df_encode[
            (df_encode['distance_to_tss'] >= -promoter_up) & 
            (df_encode['distance_to_tss'] <= promoter_down)
        ].copy()
        
        # Build list of high confidence peaks to mark them distinctively
        high_conf_peaks = pd.DataFrame()
        if 'high_confidence_site' in window_peaks.columns:
            high_conf_peaks = window_peaks[window_peaks['high_confidence_site'] == True]
        elif top_n is not None:
            if 'signal' in window_peaks.columns and window_peaks['signal'].max() > 0:
                high_conf_peaks = window_peaks.sort_values('signal', ascending=False).head(top_n)
            elif 'score' in window_peaks.columns:
                high_conf_peaks = window_peaks.sort_values('score', ascending=False).head(top_n)
        elif threshold is not None:
            if 'signal' in window_peaks.columns and window_peaks['signal'].max() > 0:
                high_conf_peaks = window_peaks[window_peaks['signal'] >= threshold]
            elif 'score' in window_peaks.columns:
                high_conf_peaks = window_peaks[window_peaks['score'] >= threshold]
                
        # Draw peaks as triangles in a multi-lane or single-lane layout
        # ChIP peaks Y range: 265 to 325. Center = 295.
        if not window_peaks.empty:
            high_conf_indices = set(high_conf_peaks.index) if not high_conf_peaks.empty else set()
            
            # Simple lane packing for overlapping peaks to make it extremely clear
            peaks_lanes = {}
            for idx, row in window_peaks.iterrows():
                rel = row['distance_to_tss']
                cx = get_x(rel)
                
                # Check for lanes
                lane = 0
                while True:
                    collision = False
                    if lane in peaks_lanes:
                        for prev_x in peaks_lanes[lane]:
                            if abs(cx - prev_x) < 14:  # minimum peak horizontal spacing
                                collision = True
                                break
                    if not collision:
                        break
                    lane += 1
                    if lane >= 3:
                        lane = 0
                        break
                
                if lane not in peaks_lanes:
                    peaks_lanes[lane] = []
                peaks_lanes[lane].append(cx)
                
                py = 273 + (lane * 16)
                
                is_high_conf = idx in high_conf_indices
                
                # Colors
                color_class = "peak-triangle-high" if is_high_conf else ""
                stroke_peak = "#ffffff" if is_high_conf else "#fda4af"
                fill_peak = "#ef4444" if is_high_conf else "#f43f5e"
                
                # SVG Path representing triangle pointing down: M cx-6,py L cx+6,py L cx,py+10 Z
                points_triangle = f"{cx-6:.1f},{py:.1f} {cx+6:.1f},{py:.1f} {cx:.1f},{py+10:.1f}"
                
                sig_val = row.get('signal', row.get('score', 0.0))
                if pd.isna(sig_val):
                    sig_val = 0.0
                
                tooltip_txt = (
                    f'<div class="tooltip-title" style="color: #fda4af;">🔥 ENCODE ChIP-Seq Peak</div>'
                    f'<div class="tooltip-row"><span class="tooltip-label">Accession:</span> <span class="tooltip-val">{row.get("file_accession", "N/A")}</span></div>'
                    f'<div class="tooltip-row"><span class="tooltip-label">TF Target:</span> <span class="tooltip-val" style="color:#ef4444; font-weight:700;">{row.get("experiment", "N/A")}</span></div>'
                    f'<div class="tooltip-row"><span class="tooltip-label">Biosample:</span> <span class="tooltip-val" style="color:#60a5fa;">{row.get("biosample", "N/A")}</span></div>'
                    f'<div class="tooltip-row"><span class="tooltip-label">Signal Intensity:</span> <span class="tooltip-val">{sig_val:.2f}</span></div>'
                    f'<div class="tooltip-row"><span class="tooltip-label">Distance to TSS:</span> <span class="tooltip-val">{rel:+,} bp</span></div>'
                    f'<div class="tooltip-row"><span class="tooltip-label">High-Confidence:</span> <span class="tooltip-val">{"✅ Yes (Conserved)" if is_high_conf else "No"}</span></div>'
                )
                
                html_out.append(f'<polygon points="{points_triangle}" class="peak-triangle {color_class}" '
                                f'style="fill: {fill_peak}; stroke: {stroke_peak};" '
                                f'onmousemove="showTooltip(event, \'{tooltip_txt}\')" onmouseout="hideTooltip()" />')
                
                # Render little signal text above high confidence peaks
                if is_high_conf and sig_val > 0:
                    html_out.append(f'<text x="{cx}" y="{py - 2}" fill="#ffffff" font-size="8px" font-weight="700" '
                                    f'text-anchor="middle" pointer-events="none">{sig_val:.1f}</text>')

    # 9. TSS Line overlay
    html_out.append(f'<line x1="{tss_x}" y1="40" x2="{tss_x}" y2="330" class="tss-line" />')
    html_out.append(f'<rect x="{tss_x - 12}" y="35" width="24" height="12" rx="3" fill="#ef4444" />')
    html_out.append(f'<text x="{tss_x}" y="44" fill="#ffffff" font-size="8px" font-weight="700" text-anchor="middle">TSS</text>')
    
    # 10. Axis Ticks & Grid (Coordinate axis line at bottom of scale Y=70)
    html_out.append(f'<line x1="{margin_l}" y1="70" x2="{margin_l + track_w}" y2="70" stroke="#475569" stroke-width="1.2" />')
    
    # Plot ticks every 250 bp or 500 bp
    tick_step = 250 if total_bp <= 3000 else 500
    # Find nice starting tick offset
    start_tick = ((-promoter_up) // tick_step) * tick_step
    while start_tick <= promoter_down:
        if start_tick >= -promoter_up:
            tx = get_x(start_tick)
            html_out.append(f'<line x1="{tx}" y1="70" x2="{tx}" y2="75" stroke="#475569" stroke-width="1.2" />')
            
            # Tick text
            label = "TSS" if start_tick == 0 else f"{start_tick:+,} bp"
            label_c = "#ef4444" if start_tick == 0 else "#64748b"
            html_out.append(f'<text x="{tx}" y="62" class="axis-label" style="fill: {label_c};" text-anchor="middle">{label}</text>')
        start_tick += tick_step
        
    # 11. Visual Legend at Bottom Y=355 to 380
    legend_y = 365
    html_out.append(f'<rect x="{margin_l}" y="{legend_y - 12}" width="{track_w}" height="24" fill="#151c27" rx="6" stroke="rgba(255,255,255,0.03)" />')
    
    # Legend Items
    # 1. TSS Line
    html_out.append(f'<line x1="{margin_l + 20}" y1="{legend_y}" x2="{margin_l + 40}" y2="{legend_y}" class="tss-line" />')
    html_out.append(f'<text x="{margin_l + 48}" y="{legend_y + 3}" class="axis-label" style="fill:#ef4444; font-weight:700;">TSS</text>')
    
    # 2. Motif Pill
    html_out.append(f'<rect x="{margin_l + 120}" y="{legend_y - 6}" width="20" height="10" rx="2" fill="#3b82f6" fill-opacity="0.6" stroke="#60a5fa" stroke-width="1" />')
    html_out.append(f'<text x="{margin_l + 146}" y="{legend_y + 3}" class="axis-label" style="fill:#60a5fa; font-weight:600;">Predicted TF Motif</text>')
    
    # 3. Peak Triangle
    html_out.append(f'<polygon points="{margin_l + 310},{legend_y - 5} {margin_l + 320},{legend_y - 5} {margin_l + 315},{legend_y + 5}" fill="#f43f5e" fill-opacity="0.7" stroke="#fda4af" stroke-width="1" />')
    html_out.append(f'<text x="{margin_l + 326}" y="{legend_y + 3}" class="axis-label" style="fill:#fda4af; font-weight:600;">ENCODE ChIP Peak</text>')
    
    # 4. Conserved Peak (High confidence)
    html_out.append(f'<polygon points="{margin_l + 480},{legend_y - 5} {margin_l + 490},{legend_y - 5} {margin_l + 485},{legend_y + 5}" fill="#ef4444" stroke="#ffffff" stroke-width="1.2" />')
    html_out.append(f'<text x="{margin_l + 496}" y="{legend_y + 3}" class="axis-label" style="fill:#ffffff; font-weight:700;">Conserved Peak</text>')
    
    # 5. Conservation Track Wave
    html_out.append(f'<rect x="{margin_l + 630}" y="{legend_y - 6}" width="20" height="10" fill="#14b8a6" fill-opacity="0.3" stroke="#14b8a6" stroke-width="1" />')
    html_out.append(f'<text x="{margin_l + 656}" y="{legend_y + 3}" class="axis-label" style="fill:#14b8a6; font-weight:600;">phastCons Score</text>')
    
    # 6. Synergy Zone
    html_out.append(f'<rect x="{margin_l + 780}" y="{legend_y - 6}" width="15" height="10" fill="#f59e0b" fill-opacity="0.25" stroke="#f59e0b" stroke-dasharray="1 1" />')
    html_out.append(f'<text x="{margin_l + 800}" y="{legend_y + 3}" class="axis-label" style="fill:#f59e0b; font-weight:700;">Synergy Zone</text>')
    
    # Close SVG
    html_out.append('</svg>')
    
    # Add HTML tooltip box
    html_out.append('<div id="svg-tooltip" class="custom-tooltip"></div>')
    
    # Close HTML wrapper
    html_out.append('</div>')
    
    # Add Hover JavaScript
    js_code = """
    <script>
        function showTooltip(evt, text) {
            var tooltip = document.getElementById('svg-tooltip');
            tooltip.innerHTML = text;
            tooltip.style.display = 'block';
            
            // Adjust positioning so it fits beautifully
            var pageX = evt.pageX;
            var pageY = evt.pageY;
            
            tooltip.style.left = (pageX + 15) + 'px';
            tooltip.style.top = (pageY - 10) + 'px';
            tooltip.style.opacity = 1;
        }
        
        function hideTooltip() {
            var tooltip = document.getElementById('svg-tooltip');
            tooltip.style.opacity = 0;
            tooltip.style.display = 'none';
        }
    </script>
    """
    
    html_out.append(js_code)
    
    return "\n".join(html_out)
