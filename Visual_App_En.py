import streamlit as st
import pandas as pd
from stmol import showmol
import py3Dmol
import os
import pathlib
import time
import logging
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import Cavity_Calculation as cc
import traj_tool
from calculate_ESP_MHP import (
    calculate_esp_for_files,
    calculate_mhp_for_files,
    detect_main_metal_charge,
)

logging.getLogger("hydrophobicity").setLevel(logging.ERROR)


# --- 1. State Management & Callback Functions ---
def clear_results():
    st.session_state.all_results = {}
    st.session_state.selected_detail_file = None


def classify_cavity_property(mean_mhp, mean_esp, hi):
    labels = []
    if mean_mhp is None:
        labels.append("MHP Undetermined")
    elif mean_mhp > 0:
        labels.append("Hydrophobic")
    elif mean_mhp < 0:
        labels.append("Hydrophilic")
    else:
        labels.append("Hydrophobic Neutral")

    if mean_esp is None:
        labels.append("ESP Undetermined")
    elif mean_esp > 0:
        labels.append("Positive Potential")
    elif mean_esp < 0:
        labels.append("Negative Potential")
    else:
        labels.append("Potential Neutral")

    if hi is not None:
        labels.append("High HI" if hi >= 0.5 else "Low HI")

    if mean_mhp is None and mean_esp is None and hi is None:
        return "Properties Not Calculated"
    return " / ".join(labels)


def format_metric(value, digits=3):
    if value is None:
        return "N/A"
    if isinstance(value, float) and np.isnan(value):
        return "N/A"
    return f"{value:.{digits}f}"


if 'cav' not in st.session_state:
    st.session_state.cav = cc.cavity()
if 'all_results' not in st.session_state:
    st.session_state.all_results = {}
if 'selected_detail_file' not in st.session_state:
    st.session_state.selected_detail_file = None

st.set_page_config(layout="wide", page_title="Supramolecular Cage Analysis System v2.3")

# --- 2. Sidebar: Configuration & Mode Switch ---
with st.sidebar:
    st.header("⚙️ Computation Console")

    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(SCRIPT_DIR, "datas", "input")
    os.environ["BABEL_DATADIR"] = data_dir
    output_dir = os.path.join(SCRIPT_DIR, "datas", "output")

    mode = st.radio(
        "Select Run Mode",
        ["Static (Single/Multi-file Analysis)", "Dynamic (Trajectory Analysis)"],
        on_change=clear_results
    )
    is_dynamic = (mode == "Dynamic (Trajectory Analysis)")

    raw_files = os.listdir(data_dir)
    files_to_process = []

    if not is_dynamic:
        pdb_list = sorted([
            f for f in raw_files
            if f.lower().endswith(('.pdb', '.mol2'))
               and not pathlib.Path(f).stem.lower().endswith('_PROBE')
        ])
        files_to_process = st.multiselect("Select PDB files for calculation", pdb_list)
    else:
        traj_list = [f for f in raw_files if not f.endswith(('.pdb', '.mol2', '.png', '.obj'))]
        selected_traj = st.selectbox("Select Trajectory File (HISTORY/XTC)", traj_list)
        if selected_traj:
            files_to_process = [selected_traj]

    ball_type = st.number_input("Ball Center Type", value=2)
    sub_time = st.number_input("Subdivision Time", value=4)

    # ================= Edge Threshold & Inheritance Control Panel =================
    if is_dynamic:
        st.markdown("---")
        dynamic_edge_threshold = st.number_input(
            "⚡ Dynamic Edge Threshold",
            value=1.0,
            step=0.1,
            format="%.2f",
            help="Controls the independent target edge length squared for mesh subdivision. Smaller values yield finer subdivisions."
        )
        enable_inheritance = st.checkbox(
            "🚀 Enable Time Inheritance Acceleration",
            value=True,
            help="[ON]: Activates loop-free lightweight engine, inheriting topological mesh tweaks from the previous frame to significantly boost frame rate; [OFF]: Breaks connection between frames, reconstructing the balloon volume entirely from scratch per frame."
        )
    else:
        with st.expander("📐 Subdivision Precision Mapping (Max Edge Length)", expanded=False):
            st.markdown("""
            Customize the **Max Edge Length (Target Val)** corresponding to different **Approximate Inner Diameters (Å)**.
            """)
            col1, col2 = st.columns(2)
            with col1:
                val_under_5 = st.number_input("Edge length for < 5 Å", value=0.05, step=0.01, format="%.2f")
                val_5_10 = st.number_input("Edge length for 5 - 10 Å", value=0.25, step=0.05, format="%.2f")
                val_10_15 = st.number_input("Edge length for 10 - 15 Å", value=0.50, step=0.05, format="%.2f")
            with col2:
                val_15_20 = st.number_input("Edge length for 15 - 20 Å", value=1.00, step=0.10, format="%.2f")
                val_over_20 = st.number_input("Edge length for > 20 Å", value=10.00, step=1.00, format="%.2f")

            edge_thresholds = {
                "under_5": val_under_5,
                "5_to_10": val_5_10,
                "10_to_15": val_10_15,
                "15_to_20": val_15_20,
                "over_20": val_over_20
            }
    # ====================================================================

    if st.button("🚀 Start Calculation"):
        st.session_state.all_results = {}
        target_list = []

        if is_dynamic:
            with st.spinner("Parsing trajectory..."):
                traj_path = os.path.join(data_dir, files_to_process[0])
                pdb_frames = traj_tool.split_trajectory_to_pdbs(
                    input_filepath=traj_path, output_directory=data_dir,
                    swap_atoms={"he": "H"}, forcefield="opls"
                )
                target_list = [pathlib.Path(p).name for p in pdb_frames]
        else:
            target_list = files_to_process

        for idx, f_name in enumerate(target_list):
            with st.spinner(f"Analyzing {f_name}..."):
                try:
                    if is_dynamic:
                        if idx == 0 and hasattr(st.session_state.cav, 'last_frame_data'):
                            st.session_state.cav.last_frame_data = None
                        elif not enable_inheritance and hasattr(st.session_state.cav, 'last_frame_data'):
                            st.session_state.cav.last_frame_data = None
                    else:
                        if hasattr(st.session_state.cav, 'last_frame_data'):
                            st.session_state.cav.last_frame_data = None

                    average_cavity_hydrophobicity = None
                    hydrophobic_index = None
                    mean_abs_esp = None
                    max_mhp = None
                    min_mhp = None
                    max_esp = None
                    min_esp = None
                    mhp_values = np.array([])
                    esp_values = np.array([])
                    property_label = "Properties Not Calculated"

                    frame_start_time = time.time()

                    if is_dynamic:
                        v_list = st.session_state.cav.Calculate_Cavity(
                            f_name, ball_type, sub_time, data_dir + "/", output_dir + "/",
                            use_dynamic_engine=enable_inheritance, dynamic_threshold=dynamic_edge_threshold
                        )
                    else:
                        st.session_state.cav.edge_thresholds = edge_thresholds
                        v_list = st.session_state.cav.Calculate_Cavity(
                            f_name, ball_type, sub_time, data_dir + "/", output_dir + "/"
                        )

                    frame_end_time = time.time()
                    frame_cost = frame_end_time - frame_start_time

                    try:
                        if v_list and v_list[0] > 10000:
                            st.warning(
                                f"⚠️ Volume of {f_name} is too large, skipping MHP/ESP grid calculation automatically.")
                        else:
                            full_pdb_path = os.path.join(data_dir, f_name)
                            stem = pathlib.Path(f_name).stem

                            probe_path = os.path.join(output_dir, f"{stem}_cavity.pdb")

                            if not os.path.exists(probe_path):
                                st.warning(
                                    f"⚠️ {f_name} is missing the cavity file, cannot calculate MHP/ESP: {probe_path}")
                            else:
                                try:
                                    metal, metal_charge = detect_main_metal_charge(full_pdb_path)
                                    esp_result = calculate_esp_for_files(
                                        full_pdb_path,
                                        probe_path,
                                        method="eem",
                                        metal_name=metal,
                                        metal_charge=metal_charge,
                                    )
                                    esp_values = np.asarray(esp_result["esp"], dtype=float)
                                    mean_abs_esp = float(np.mean(esp_values))
                                    max_esp = float(np.max(esp_values))
                                    min_esp = float(np.min(esp_values))
                                except Exception as esp_e:
                                    st.warning(f"⚠️ ESP calculation failed for {f_name}: {esp_e}")

                                try:
                                    mhp_result = calculate_mhp_for_files(
                                        full_pdb_path,
                                        probe_path,
                                        method="Ghose",
                                        distance_function="Fauchere",
                                    )
                                    mhp_values = np.asarray(mhp_result["mhp"], dtype=float)
                                    average_cavity_hydrophobicity = float(np.mean(mhp_values))
                                    max_mhp = float(np.max(mhp_values))
                                    min_mhp = float(np.min(mhp_values))

                                    mlp_pos = mhp_values[mhp_values > 0]
                                    mlp_neg = mhp_values[mhp_values < 0]
                                    denom = float(np.sum(mlp_pos) - np.sum(mlp_neg))
                                    hydrophobic_index = float(np.sum(mlp_pos) / denom) if denom != 0 else 0
                                except Exception as mhp_e:
                                    st.warning(
                                        f"⚠️ MHP calculation failed for {f_name}, ESP results retained: {mhp_e}")
                    except Exception as chem_e:
                        st.warning(f"⚠️ Physicochemical property analysis failed for {f_name}: {chem_e}")

                    property_label = classify_cavity_property(
                        average_cavity_hydrophobicity,
                        mean_abs_esp,
                        hydrophobic_index,
                    )

                    if v_list:
                        st.session_state.all_results[f_name] = {
                            'vol': v_list[0],
                            'cost_time': frame_cost,
                            'win_info': st.session_state.cav.last_window_info,
                            'rebek': st.session_state.cav.cached_rebek_vol,
                            'pdb_data': open(os.path.join(data_dir, f_name), 'r').read(),
                            'mhp': average_cavity_hydrophobicity,
                            'max_mhp': max_mhp,
                            'min_mhp': min_mhp,
                            'hi': hydrophobic_index,
                            'mean_esp': mean_abs_esp,
                            'max_esp': max_esp,
                            'min_esp': min_esp,
                            'mhp_grid': mhp_values.tolist() if mhp_values.size > 0 else None,
                            'esp_grid': esp_values.tolist() if esp_values.size > 0 else None,
                            'property_label': property_label
                        }
                except Exception as e:
                    st.error(f"❌ Critical error occurred while analyzing {f_name}: {e}")

# --- 3. Data Processing & Summary Charts ---
if st.session_state.all_results:
    summary_list = []
    abs_errors = []
    window_data_dict = {}

    for name, data in st.session_state.all_results.items():
        calc_v = data.get('vol', 0)
        rebek_v = data.get('rebek')

        row = {
            "File Name": name,
            "Calculated Vol (Å³)": round(calc_v, 2) if calc_v else "N/A",
        }

        if is_dynamic:
            row["Cost Time (s)"] = round(data.get('cost_time', 0), 3)

        diff_str = "N/A"
        if rebek_v and rebek_v != 0 and calc_v:
            err = ((calc_v - rebek_v) / rebek_v) * 100
            diff_str = f"{err:+.2f}%"
            abs_errors.append(abs(err))

        row["Theoretical Vol (Rebek)"] = round(rebek_v, 2) if rebek_v else "N/A"
        row["Relative Error"] = diff_str

        if data.get('win_info'):
            win_data = data['win_info']
            row["Max Window"] = round(
                max([w['diameter'] for w in win_data.get('windows', [])]) if win_data.get('windows') else 0, 2)
            row["Window Count"] = win_data.get('window_count', 0)
            cage_name = name.split('.')[0]
            window_data_dict[cage_name] = [w['diameter'] for w in win_data.get('windows', [])]

        row["Avg Hydrophobicity (MHP)"] = round(data.get('mhp'), 4) if data.get('mhp') is not None else "N/A"
        row["Min MHP"] = round(data.get('min_mhp'), 4) if data.get('min_mhp') is not None else "N/A"
        row["Max MHP"] = round(data.get('max_mhp'), 4) if data.get('max_mhp') is not None else "N/A"
        row["Hydrophobic Index (HI)"] = round(data.get('hi'), 3) if data.get('hi') is not None else "N/A"
        row["Avg ESP"] = round(data.get('mean_esp'), 2) if data.get('mean_esp') is not None else "N/A"
        row["Min ESP"] = round(data.get('min_esp'), 2) if data.get('min_esp') is not None else "N/A"
        row["Max ESP"] = round(data.get('max_esp'), 2) if data.get('max_esp') is not None else "N/A"
        row["Determined Properties"] = data.get('property_label', "Properties Not Calculated")

        summary_list.append(row)

    df_summary = pd.DataFrame(summary_list)

    if is_dynamic:
        df_summary['frame_idx'] = df_summary['File Name'].str.extract(r'(\d+)').astype(float)
        df_summary = df_summary.sort_values('frame_idx')

        st.subheader("📈 Trajectory Volume Evolution & Cost Time Monitoring (Dual-Axis)")

        fig_dual = make_subplots(specs=[[{"secondary_y": True}]])

        fig_dual.add_trace(go.Bar(
            x=df_summary["File Name"],
            y=df_summary["Calculated Vol (Å³)"],
            name="Cavity Vol (Å³)",
            marker_color='rgb(254, 179, 174)',
            text=df_summary["Calculated Vol (Å³)"].round(2),
            textposition='outside',
            textfont=dict(size=14, color='#000000', family="Arial Black")
        ), secondary_y=False)

        fig_dual.add_trace(go.Scatter(
            x=df_summary["File Name"],
            y=df_summary["Cost Time (s)"],
            name="Frame Time (s)",
            mode='lines+markers+text',
            line=dict(color='rgb(21, 151, 165)', width=3.5),
            marker=dict(size=9, symbol="circle", color='rgb(15, 120, 130)'),
            text=df_summary["Cost Time (s)"].round(3),
            textposition='bottom right',
            textfont=dict(size=14, color='#111111', family="Arial Black")
        ), secondary_y=True)

        fig_dual.update_layout(
            xaxis=dict(
                title=dict(text="<b>Trajectory Frame File</b>", font=dict(size=14, color="#111111")),
                tickfont=dict(size=13, color="#111111")
            ),
            yaxis=dict(
                title=dict(text="<b>Cavity Volume (Å³)</b>", font=dict(size=14, color="rgb(254, 179, 174)")),
                tickfont=dict(size=13, color="rgb(254, 179, 174)"),
                gridcolor="rgba(200, 200, 200, 0.25)"
            ),
            yaxis2=dict(
                title=dict(text="<b>Calculation Time (Seconds)</b>", font=dict(size=14, color="rgb(21, 151, 165)")),
                tickfont=dict(size=13, color="rgb(21, 151, 165)"),
                anchor="x",
                overlaying="y",
                side="right"
            ),
            legend=dict(
                x=1.02,
                y=1.0,
                xanchor="left",
                yanchor="top",
                orientation="v",
                font=dict(size=13, color="#111111")
            ),
            margin=dict(l=60, r=130, t=50, b=50),
            height=540,
            paper_bgcolor="white"
        )

        st.plotly_chart(fig_dual, use_container_width=True)
        df_summary = df_summary.drop(columns=['frame_idx'], errors='ignore')
    else:
        if abs_errors:
            mape = np.mean(abs_errors)
            st.info(
                f"📋 **Statistical Summary**: Processed {len(summary_list)} files | **Mean Absolute Percentage Error (MAPE): {mape:.2f}%**")
        if window_data_dict:
            with open(os.path.join(output_dir, "window_data.json"), 'w', encoding='utf-8') as f_json:
                import json

                json.dump(window_data_dict, f_json, indent=4)

    st.subheader("📋 Data Overview")
    event = st.dataframe(df_summary, width="stretch", hide_index=True,
                         on_select="rerun", selection_mode="single-row")

    if event.selection.rows:
        st.session_state.selected_detail_file = df_summary.iloc[event.selection.rows[0]]["File Name"]

# --- 4. Detail Linkage Area ---
if st.session_state.selected_detail_file:
    f_name = st.session_state.selected_detail_file
    res = st.session_state.all_results.get(f_name)
    if res:
        st.divider()
        main_left, main_right = st.columns([1.2, 1])

        val_hi = res.get('hi')
        val_mhp = res.get('mhp')
        val_max_mhp = res.get('max_mhp')
        val_min_mhp = res.get('min_mhp')
        val_esp = res.get('mean_esp')
        val_max_esp = res.get('max_esp')
        val_min_esp = res.get('min_esp')
        val_property = res.get('property_label', "Properties Not Calculated")
        val_vol = res.get('vol', 0)

        # 提取窗口数量用于雷达图
        val_win_count = 0
        if res.get('win_info'):
            val_win_count = res['win_info'].get('window_count', 0)

        with main_left:
            st.subheader(f"🌐 3D Model Render Observation: {f_name}")

            viz_mode = st.radio("Select Analysis Dimension",
                                ["ESP (Electrostatic Potential)", "MHP (Hydrophobicity)"], horizontal=True)
            prop_key = "ESP" if "ESP" in viz_mode else "MHP"
            current_prop_val = val_esp if prop_key == "ESP" else val_mhp
            is_prop_available = current_prop_val is not None

            view = py3Dmol.view(width='100%', height=400)
            view.addModel(res['pdb_data'], 'pdb')
            view.setStyle({'stick': {'colorscheme': 'greenCarbon'}, 'sphere': {'scale': 0.3}})

            if is_prop_available:
                surf_color = '#ff4d4d' if current_prop_val > 0 else '#4d4dff'
                view.addSurface(py3Dmol.VDW, {'opacity': 0.6, 'color': surf_color})
                view.zoomTo()
                showmol(view, height=400)

                l_min, l_max, colors = ("Negative", "Positive", "blue, white, red") if prop_key == "ESP" else (
                    "Hydrophilic", "Hydrophobic", "#3333ff, white, #ff3333")
                st.markdown(f"""
                            <div style="display: flex; align-items: center; justify-content: center; margin-top: -10px;">
                                <span style="font-size: 0.8rem; width: 80px; text-align: right; margin-right: 15px;">{l_min}</span>
                                <div style="width: 100%; max-width: 320px; height: 12px; background: linear-gradient(to right, {colors}); border-radius: 10px; border: 1px solid #ddd;"></div>
                                <span style="font-size: 0.8rem; width: 80px; margin-left: 15px;">{l_max}</span>
                            </div>
                        """, unsafe_allow_html=True)
            else:
                view.zoomTo()
                showmol(view, height=400)
                st.warning(f"⚠️ Property {prop_key} not calculated, cannot display mesh surface.")

            st.write("")
            st.divider()

            if prop_key == "MHP":
                st.markdown("**💧 Hydrophobicity (MHP) Metrics**")
                c1, c2, c3, c4 = st.columns(4)
                c1.metric("HI (Hydrophobic Index)", format_metric(val_hi, 3))
                c2.metric("MHP (Average)", format_metric(val_mhp, 3))
                c3.metric("Min MHP", format_metric(val_min_mhp, 3))
                c4.metric("Max MHP", format_metric(val_max_mhp, 3))
            else:
                st.markdown("**⚡ Electrostatic Potential (ESP) Metrics**")
                c1, c2, c3 = st.columns(3)
                c1.metric("ESP (Average)", format_metric(val_esp, 2))
                c2.metric("Min ESP", format_metric(val_min_esp, 2))
                c3.metric("Max ESP", format_metric(val_max_esp, 2))

            st.write("")
            st.info(f"Final Determined Properties: {val_property}")

            st.write(f"**📈 {prop_key} Data Distribution Histogram**")
            hist_data = res.get('esp_grid') if prop_key == "ESP" else res.get('mhp_grid')

            if hist_data is not None:
                fig_hist = px.histogram(
                    hist_data, nbins=35,
                    color_discrete_sequence=['#636EFA'] if "ESP" in viz_mode else ['#EF553B'], opacity=0.75
                )
                fig_hist.update_layout(height=180, margin=dict(l=0, r=0, t=0, b=0), showlegend=False,
                                       xaxis_title=f"{prop_key} Value", yaxis_title="Count")
                st.plotly_chart(fig_hist, use_container_width=True)

        with main_right:
            if res.get('win_info'):
                with st.expander("📊 View Window Details", expanded=True):
                    st.dataframe(
                        pd.DataFrame(res['win_info']['windows'])[['id', 'diameter', 'vertex_count']],
                        hide_index=True, width="stretch"
                    )

            if is_dynamic:
                st.write("**⏱️ Single Frame Efficiency Metrics**")
                st.metric("Loop-free engine cost", f"{res.get('cost_time', 0):.3f} s / frame")

            # 🔴 核心修改点：更新雷达图的维度
            st.write("**🕸️ Cavity Comprehensive Properties Radar Chart**")

            # 对各变量进行 0~1 的归一化，以便在同一个雷达图上显示
            r_vol = val_vol / 600
            r_esp = abs(val_esp) / 50 if val_esp is not None else 0
            r_mhp = (val_mhp + 0.5) / 1.0 if val_mhp is not None else 0

            # 设定 10 个窗口为满分 (1.0)，使用 min 函数防止溢出图表边界
            r_win = min(val_win_count / 10.0, 1.0)

            fig_radar = go.Figure(data=go.Scatterpolar(
                r=[r_vol, r_esp, r_mhp, r_win],
                theta=['Volume', 'Avg ESP', 'Avg MHP', 'Window Count'],
                fill='toself', line=dict(color='#FF4B4B')
            ))
            fig_radar.update_layout(
                polar=dict(radialaxis=dict(visible=True, range=[0, 1])),
                height=320, margin=dict(l=40, r=40, t=40, b=40)
            )
            st.plotly_chart(fig_radar, use_container_width=True)