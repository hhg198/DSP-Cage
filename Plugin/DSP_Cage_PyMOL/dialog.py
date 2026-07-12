"""Qt interface for the DSP-Cage PyMOL plugin."""

from pathlib import Path

from pymol import cmd
from pymol.Qt import QtWidgets

from .backend import run_calculation


class DSPCageDialog(QtWidgets.QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.setWindowTitle("DSP-Cage")
        self.resize(560, 460)

        self.object_box = QtWidgets.QComboBox()

        self.refresh_button = QtWidgets.QPushButton("刷新")
        self.refresh_button.clicked.connect(self.refresh_objects)

        self.output_edit = QtWidgets.QLineEdit(
            str(Path.home() / "DSP_Cage_results")
        )

        self.output_button = QtWidgets.QPushButton("选择")
        self.output_button.clicked.connect(self.select_output_directory)

        self.center_type_box = QtWidgets.QComboBox()
        self.center_type_box.addItem("Type 1", 1)
        self.center_type_box.addItem("Type 2", 2)
        self.center_type_box.addItem("Type 3", 3)
        self.center_type_box.setCurrentIndex(1)

        self.subdivision_spin = QtWidgets.QSpinBox()
        self.subdivision_spin.setRange(1, 8)
        self.subdivision_spin.setValue(4)

        self.run_button = QtWidgets.QPushButton("Calculate Cavity")
        self.run_button.clicked.connect(self.run_analysis)

        self.result_text = QtWidgets.QPlainTextEdit()
        self.result_text.setReadOnly(True)

        self.close_button = QtWidgets.QPushButton("关闭")
        self.close_button.clicked.connect(self.hide)

        self.build_layout()
        self.refresh_objects()

    def build_layout(self):
        object_widget = QtWidgets.QWidget()
        object_layout = QtWidgets.QHBoxLayout(object_widget)
        object_layout.setContentsMargins(0, 0, 0, 0)
        object_layout.addWidget(self.object_box)
        object_layout.addWidget(self.refresh_button)

        output_widget = QtWidgets.QWidget()
        output_layout = QtWidgets.QHBoxLayout(output_widget)
        output_layout.setContentsMargins(0, 0, 0, 0)
        output_layout.addWidget(self.output_edit)
        output_layout.addWidget(self.output_button)

        form_layout = QtWidgets.QFormLayout()
        form_layout.addRow("PyMOL object:", object_widget)
        form_layout.addRow("Output directory:", output_widget)
        form_layout.addRow("Center type:", self.center_type_box)
        form_layout.addRow("Subdivision:", self.subdivision_spin)

        button_layout = QtWidgets.QHBoxLayout()
        button_layout.addWidget(self.run_button)
        button_layout.addWidget(self.close_button)

        main_layout = QtWidgets.QVBoxLayout(self)
        main_layout.addLayout(form_layout)
        main_layout.addLayout(button_layout)
        main_layout.addWidget(QtWidgets.QLabel("Calculation results:"))
        main_layout.addWidget(self.result_text)

    def refresh_objects(self):
        current_name = self.object_box.currentText()

        self.object_box.clear()

        object_names = cmd.get_object_list()
        self.object_box.addItems(object_names)

        if current_name in object_names:
            self.object_box.setCurrentText(current_name)

    def select_output_directory(self):
        directory = QtWidgets.QFileDialog.getExistingDirectory(
            self,
            "Select output directory",
            self.output_edit.text(),
        )

        if directory:
            self.output_edit.setText(directory)

    def run_analysis(self):
        object_name = self.object_box.currentText().strip()

        if not object_name:
            QtWidgets.QMessageBox.warning(
                self,
                "DSP-Cage",
                "PyMOL 中没有可用于计算的分子对象。",
            )
            return

        self.run_button.setEnabled(False)
        self.run_button.setText("Calculating...")
        self.result_text.setPlainText(
            f"正在计算：{object_name}\n"
            "计算过程中请勿关闭 PyMOL。"
        )

        QtWidgets.QApplication.processEvents()

        try:
            result = run_calculation(
                selection=object_name,
                output_dir=self.output_edit.text().strip(),
                center_type=self.center_type_box.currentData(),
                subdivision=self.subdivision_spin.value(),
            )

            window_info = result.get("window_info") or {}
            windows = window_info.get("windows", [])

            result_lines = [
                "DSP-Cage calculation completed",
                "",
                f"Object: {object_name}",
                f"Atom count: {result['atom_count']}",
                f"Cavity volume: {result['volume']:.3f} Å³",
                f"Window count: {window_info.get('window_count', 0)}",
            ]

            rebek_volume = result.get("rebek_volume")

            if rebek_volume is not None:
                result_lines.append(
                    f"Rebek volume: {float(rebek_volume):.3f} Å³"
                )

            if windows:
                result_lines.append("")
                result_lines.append("Window diameters:")

                for window in windows:
                    result_lines.append(
                        "  Window "
                        f"{window.get('id', '-')}: "
                        f"{window.get('diameter', 0):.3f} Å"
                    )

            result_lines.extend(
                [
                    "",
                    f"Input file: {result['input_file']}",
                    f"Cavity file: {result.get('cavity_file') or 'Not generated'}",
                    f"Output directory: {result['output_dir']}",
                ]
            )

            self.result_text.setPlainText("\n".join(result_lines))

        except Exception as exc:
            self.result_text.setPlainText(
                f"Calculation failed\n\n"
                f"{type(exc).__name__}: {exc}"
            )

            QtWidgets.QMessageBox.critical(
                self,
                "DSP-Cage calculation failed",
                str(exc),
            )

        finally:
            self.run_button.setEnabled(True)
            self.run_button.setText("Calculate Cavity")