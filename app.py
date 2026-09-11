import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from io import BytesIO

# RDKit
from rdkit import Chem
from rdkit.Chem import (
    Draw,
    AllChem,
    Descriptors,
    Crippen,
    Lipinski,
    rdMolDescriptors
)

# 3D Visualization
import py3Dmol
import streamlit.components.v1 as components

# PDF
from reportlab.lib.pagesizes import A4
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY
from reportlab.lib.colors import HexColor
from reportlab.platypus import (
    SimpleDocTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
    PageBreak,
    Image
)
from reportlab.lib.units import inch


# ============================================================
# PAGE CONFIGURATION
# ============================================================

st.set_page_config(
    page_title="Cheminformatics Virtual Laboratory",
    page_icon="🧪",
    layout="wide"
)


# ============================================================
# CUSTOM CSS
# ============================================================

st.markdown("""
<style>

.main-title {
    text-align: center;
    font-size: 40px;
    font-weight: bold;
}

.section-title {
    font-size: 25px;
    font-weight: bold;
}

</style>
""", unsafe_allow_html=True)


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def get_molecule(smiles):
    """Convert SMILES into an RDKit molecule."""

    if not smiles or not smiles.strip():
        return None

    try:
        mol = Chem.MolFromSmiles(smiles.strip())
        return mol

    except Exception:
        return None


# ============================================================
# CALCULATE COMPLETE MOLECULAR PROPERTIES
# ============================================================

def calculate_properties(mol):

    # Calculate bond information
    bond_info = calculate_bond_information(mol)

    molecular_formula = rdMolDescriptors.CalcMolFormula(mol)

    properties = {

        "Molecular Formula": molecular_formula,

        "Molecular Weight": round(
            Descriptors.MolWt(mol),
            3
        ),

        "Exact Molecular Weight": round(
            Descriptors.ExactMolWt(mol),
            4
        ),

        "LogP": round(
            Crippen.MolLogP(mol),
            3
        ),

        "TPSA": round(
            rdMolDescriptors.CalcTPSA(mol),
            3
        ),

        "Hydrogen Bond Donors": Lipinski.NumHDonors(mol),

        "Hydrogen Bond Acceptors": Lipinski.NumHAcceptors(mol),

        "Rotatable Bonds": Lipinski.NumRotatableBonds(mol),

        "Ring Count": Lipinski.RingCount(mol),

        "Aromatic Rings": rdMolDescriptors.CalcNumAromaticRings(mol),

        "Aliphatic Rings": rdMolDescriptors.CalcNumAliphaticRings(mol),

        "Number of Atoms": mol.GetNumAtoms(),

        "Heavy Atoms": mol.GetNumHeavyAtoms(),

        # Bond information
        "Heavy Atom Bonds": bond_info["Heavy Atom Bonds"],

        "Sigma (σ) Bonds": bond_info["Sigma Bonds"],

        "Pi (π) Bonds": bond_info["Pi Bonds"],

        "Total Bonds (σ + π)": bond_info["Total Bonds"],

        "Fraction Csp3": round(
            rdMolDescriptors.CalcFractionCSP3(mol),
            3
        ),

        "Molar Refractivity": round(
            Crippen.MolMR(mol),
            3
        )
    }

    return properties

# ============================================================
# GENERATE 2D IMAGE
# ============================================================

def generate_2d_image(mol):

    image = Draw.MolToImage(
        mol,
        size=(700, 500)
    )

    return image


# ============================================================
# GENERATE 3D MOLECULE
# ============================================================

def generate_3d_molecule(mol):

    try:

        mol3d = Chem.AddHs(mol)

        params = AllChem.ETKDGv3()

        params.randomSeed = 42

        status = AllChem.EmbedMolecule(
            mol3d,
            params
        )

        if status != 0:
            return None

        try:

            AllChem.MMFFOptimizeMolecule(
                mol3d
            )

        except Exception:

            try:

                AllChem.UFFOptimizeMolecule(
                    mol3d
                )

            except Exception:
                pass

        return mol3d

    except Exception:

        return None


# ============================================================
# DISPLAY 3D STRUCTURE
# ============================================================

def display_3d_structure(mol3d, style, show_surface=False):

    mol_block = Chem.MolToMolBlock(mol3d)

    view = py3Dmol.view(
        width=900,
        height=600
    )

    view.addModel(
        mol_block,
        "mol"
    )

    # --------------------------------------------------------
    # BALL AND STICK
    # --------------------------------------------------------

    if style == "Ball and Stick":

        view.setStyle({

            "stick": {
                "radius": 0.15
            },

            "sphere": {
                "scale": 0.30
            }

        })

    # --------------------------------------------------------
    # STICK
    # --------------------------------------------------------

    elif style == "Stick":

        view.setStyle({

            "stick": {
                "radius": 0.20
            }

        })

    # --------------------------------------------------------
    # SPACE FILLING
    # --------------------------------------------------------

    elif style == "Space Filling":

        view.setStyle({

            "sphere": {
                "scale": 1.0
            }

        })

    # --------------------------------------------------------
    # WIREFRAME
    # --------------------------------------------------------

    elif style == "Wireframe":

        view.setStyle({

            "line": {
                "linewidth": 2
            }

        })

    # --------------------------------------------------------
    # CARTOON STYLE
    # --------------------------------------------------------

    elif style == "Cartoon":

        # Cartoon is mainly useful for biomolecules,
        # but this provides a molecular representation.

        view.setStyle({

            "stick": {},

            "sphere": {
                "scale": 0.20
            }

        })

    # --------------------------------------------------------
    # OPTIONAL SURFACE
    # --------------------------------------------------------

    if show_surface:

        view.addSurface(
            py3Dmol.VDW,
            {
                "opacity": 0.7
            }
        )

    view.setBackgroundColor("white")

    view.zoomTo()

    html = view._make_html()

    components.html(
        html,
        height=620
    )


# ============================================================
# CREATE MOLECULAR PDF REPORT
# ============================================================

def generate_molecule_pdf(
    input_smiles,
    mol,
    properties
):

    buffer = BytesIO()

    doc = SimpleDocTemplate(
        buffer,
        pagesize=A4,
        leftMargin=45,
        rightMargin=45,
        topMargin=45,
        bottomMargin=45
    )

    styles = getSampleStyleSheet()

    title_style = ParagraphStyle(
        "TitleStyle",
        parent=styles["Title"],
        fontSize=20,
        alignment=TA_CENTER,
        textColor=HexColor("#17365D")
    )

    heading_style = ParagraphStyle(
        "HeadingStyle",
        parent=styles["Heading1"],
        fontSize=14,
        textColor=HexColor("#17365D")
    )

    body_style = ParagraphStyle(
        "BodyStyle",
        parent=styles["BodyText"],
        fontSize=10,
        leading=14,
        alignment=TA_JUSTIFY
    )

    story = []

    story.append(
        Paragraph(
            "MOLECULAR PROPERTIES REPORT",
            title_style
        )
    )

    story.append(
        Spacer(1, 0.3 * inch)
    )

    # --------------------------------------------------------
    # BASIC INFORMATION
    # --------------------------------------------------------

    story.append(
        Paragraph(
            "1. Molecular Information",
            heading_style
        )
    )

    canonical_smiles = Chem.MolToSmiles(mol)

    basic_data = [

        ["Property", "Value"],

        ["Input SMILES", input_smiles],

        ["Canonical SMILES", canonical_smiles],

        [
            "Molecular Formula",
            properties["Molecular Formula"]
        ]
    ]

    basic_table = Table(
        basic_data,
        colWidths=[2 * inch, 4 * inch]
    )

    basic_table.setStyle(
        TableStyle([

            ("BACKGROUND", (0, 0), (-1, 0), HexColor("#17365D")),

            ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),

            ("GRID", (0, 0), (-1, -1), 0.5, colors.grey),

            ("VALIGN", (0, 0), (-1, -1), "TOP"),

            ("PADDING", (0, 0), (-1, -1), 7)

        ])
    )

    story.append(basic_table)

    story.append(
        Spacer(1, 0.3 * inch)
    )

    # --------------------------------------------------------
    # PROPERTY TABLE
    # --------------------------------------------------------

    story.append(
        Paragraph(
            "2. Molecular Properties",
            heading_style
        )
    )

    property_data = [
        ["Property", "Value"]
    ]

    for key, value in properties.items():

        property_data.append(
            [key, str(value)]
        )

    property_table = Table(
        property_data,
        colWidths=[3 * inch, 3 * inch],
        repeatRows=1
    )

    property_table.setStyle(
        TableStyle([

            ("BACKGROUND", (0, 0), (-1, 0), HexColor("#2F5597")),

            ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),

            ("GRID", (0, 0), (-1, -1), 0.5, colors.grey),

            ("VALIGN", (0, 0), (-1, -1), "TOP"),

            ("PADDING", (0, 0), (-1, -1), 6)

        ])
    )

    story.append(property_table)

    story.append(
        Spacer(1, 0.3 * inch)
    )

    story.append(
        Paragraph(
            "3. Interpretation",
            heading_style
        )
    )

    interpretation = """
    This report was generated automatically using the
    Cheminformatics Virtual Laboratory. The molecular properties
    were calculated from the chemical structure represented by
    the submitted SMILES notation using the RDKit cheminformatics
    toolkit.
    """

    story.append(
        Paragraph(
            interpretation,
            body_style
        )
    )

    # --------------------------------------------------------
    # BUILD PDF
    # --------------------------------------------------------

    doc.build(story)

    pdf_data = buffer.getvalue()

    buffer.close()

    return pdf_data


# ============================================================
# GENERATE GENERAL FINAL PROJECT REPORT
# ============================================================

def generate_final_report():

    buffer = BytesIO()

    doc = SimpleDocTemplate(
        buffer,
        pagesize=A4,
        leftMargin=50,
        rightMargin=50,
        topMargin=55,
        bottomMargin=50
    )

    styles = getSampleStyleSheet()

    title_style = ParagraphStyle(
        "ProjectTitle",
        parent=styles["Title"],
        fontSize=22,
        leading=28,
        alignment=TA_CENTER,
        textColor=HexColor("#17365D")
    )

    subtitle_style = ParagraphStyle(
        "ProjectSubtitle",
        parent=styles["Normal"],
        fontSize=13,
        leading=18,
        alignment=TA_CENTER
    )

    heading_style = ParagraphStyle(
        "ProjectHeading",
        parent=styles["Heading1"],
        fontSize=16,
        leading=20,
        textColor=HexColor("#17365D")
    )

    body_style = ParagraphStyle(
        "ProjectBody",
        parent=styles["BodyText"],
        fontSize=10,
        leading=15,
        alignment=TA_JUSTIFY
    )

    story = []

    # COVER PAGE

    story.append(
        Spacer(1, 1 * inch)
    )

    story.append(
        Paragraph(
            "FINAL PROJECT REPORT",
            subtitle_style
        )
    )

    story.append(
        Spacer(1, 0.3 * inch)
    )

    story.append(
        Paragraph(
            "Cheminformatics Virtual Laboratory",
            title_style
        )
    )

    story.append(
        Spacer(1, 0.3 * inch)
    )

    story.append(
        Paragraph(
            "An Interactive Platform for Molecular Representation, "
            "Visualization and Molecular Property Analysis",
            subtitle_style
        )
    )

    story.append(
        PageBreak()
    )

    sections = [

        (
            "1. Abstract",

            "The Cheminformatics Virtual Laboratory was developed "
            "as an interactive web-based platform for teaching "
            "molecular representation, visualization and molecular "
            "property calculation."
        ),

        (
            "2. Objectives",

            "The objectives include molecular visualization, SMILES "
            "representation, descriptor calculation, structure-property "
            "analysis and student assessment."
        ),

        (
            "3. Technologies Used",

            "The project uses Python, Streamlit, RDKit, Pandas, "
            "Matplotlib, py3Dmol and ReportLab."
        ),

        (
            "4. Functional Modules",

            "The application contains Home, Theory, Molecular "
            "Visualization, Molecular Descriptor Calculation, "
            "Structure–Property Analysis, Assessment and Final Report modules."
        ),

        (
            "5. Conclusion",

            "The application provides an interactive virtual laboratory "
            "environment for students to explore molecular structures "
            "and their physicochemical properties."
        )
    ]

    for heading, content in sections:

        story.append(
            Paragraph(
                heading,
                heading_style
            )
        )

        story.append(
            Spacer(1, 0.1 * inch)
        )

        story.append(
            Paragraph(
                content,
                body_style
            )
        )

        story.append(
            Spacer(1, 0.2 * inch)
        )

    doc.build(story)

    pdf_data = buffer.getvalue()

    buffer.close()

    return pdf_data


# ============================================================
# SIDEBAR
# ============================================================

st.sidebar.title("🧪 Virtual Lab")

page = st.sidebar.radio(

    "Navigation",

    [

        "🏠 Home",

        "📖 Theory",

        "🧬 Molecular Visualization",

        "📊 Molecular Properties Report",

        "📈 Structure–Property Analysis",

        "📝 Assessment",

        "📄 Final Project Report"

    ]
)

st.sidebar.divider()

st.sidebar.markdown("""
### 🎓 Learning Objectives

Students can learn:

- SMILES notation
- Molecular representation
- 2D visualization
- 3D visualization
- Molecular descriptors
- Structure-property relationships
""")


# ============================================================
# HOME PAGE
# ============================================================

if page == "🏠 Home":

    st.markdown(
        '<div class="main-title">'
        '🧪 Cheminformatics Virtual Laboratory'
        '</div>',
        unsafe_allow_html=True
    )

    st.markdown("## Welcome to the Virtual Laboratory")

    st.write("""
    This interactive virtual laboratory helps students understand
    important concepts in cheminformatics and computational chemistry.
    """)

    st.markdown("### 🧪 What can students do?")

    st.markdown("""
    ✅ Enter any SMILES notation

    ✅ Validate molecular structures

    ✅ Generate 2D structures

    ✅ Generate interactive 3D structures

    ✅ Select different molecular visualization styles

    ✅ Calculate complete molecular properties

    ✅ Download molecular reports

    ✅ Compare multiple molecules

    ✅ Complete an assessment
    """)

    st.markdown("### Example Molecules")

    example_df = pd.DataFrame({

        "Molecule": [

            "Ethanol",

            "Benzene",

            "Phenol",

            "Acetic Acid",

            "Aspirin",

            "Caffeine"
        ],

        "SMILES": [

            "CCO",

            "c1ccccc1",

            "Oc1ccccc1",

            "CC(=O)O",

            "CC(=O)Oc1ccccc1C(=O)O",

            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
        ]
    })

    st.dataframe(
        example_df,
        use_container_width=True,
        hide_index=True
    )


# ============================================================
# THEORY PAGE
# ============================================================

elif page == "📖 Theory":

    st.title("📖 Introduction to Cheminformatics")

    st.markdown("""
## What is Cheminformatics?

Cheminformatics combines:

- Chemistry ⚗️
- Computer Science 💻
- Data Analysis 📊

It uses computational techniques to store, process,
visualize and analyze chemical information.
""")

    st.markdown("""
## SMILES

**SMILES** means:

### Simplified Molecular Input Line Entry System

It represents a molecular structure using a text string.

Examples:

| Molecule | SMILES |
|---|---|
| Ethanol | CCO |
| Benzene | c1ccccc1 |
| Acetic Acid | CC(=O)O |
""")

    st.markdown("""
## Molecular Descriptors

Molecular descriptors provide numerical information
about a chemical structure.

Examples include:

- Molecular Weight
- LogP
- TPSA
- Hydrogen Bond Donors
- Hydrogen Bond Acceptors
- Rotatable Bonds
""")

    st.markdown("""
## 3D Molecular Models

Different visualization models help students understand
the molecular geometry.

### 🟢 Ball and Stick

Shows atoms as spheres and bonds as sticks.

### 🟡 Stick Model

Emphasizes chemical bonding.

### 🔵 Space Filling

Shows the approximate molecular volume.

### ⚫ Wireframe

Provides a simplified representation of the molecule.
""")


# ============================================================
# MOLECULAR VISUALIZATION PAGE
# ============================================================

elif page == "🧬 Molecular Visualization":

    st.title("🧬 Molecular Visualization Laboratory")

    examples = {

        "Ethanol": "CCO",

        "Benzene": "c1ccccc1",

        "Phenol": "Oc1ccccc1",

        "Acetic Acid": "CC(=O)O",

        "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",

        "Caffeine":
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    }

    col1, col2 = st.columns([2, 1])

    with col1:

        smiles = st.text_input(
            "Enter SMILES",
            value="CCO",
            key="visualization_smiles"
        )

    with col2:

        selected = st.selectbox(
            "Choose Example Molecule",
            list(examples.keys())
        )

        if st.button("Load Example"):

            smiles = examples[selected]

            st.session_state.visualization_smiles = smiles

            st.rerun()

    mol = get_molecule(smiles)

    if mol is None:

        st.error("❌ Invalid SMILES notation.")

    else:

        st.success("✅ Valid Molecular Structure!")

        col1, col2 = st.columns([1.2, 1])

        with col1:

            st.subheader("2D Molecular Structure")

            image = generate_2d_image(mol)

            st.image(image)

        with col2:

            st.subheader("Molecular Information")

            st.write(
                "**Input SMILES:**",
                smiles
            )

            st.write(
                "**Canonical SMILES:**",
                Chem.MolToSmiles(mol)
            )

            st.write(
                "**Formula:**",
                rdMolDescriptors.CalcMolFormula(mol)
            )

            st.write(
                "**Atoms:**",
                mol.GetNumAtoms()
            )

            st.write(
                "**Bonds:**",
                mol.GetNumBonds()
            )

        st.divider()

        st.subheader("🌐 Interactive 3D Molecular Structure")

        visualization_style = st.selectbox(

            "Select Molecular Model",

            [

                "Ball and Stick",

                "Stick",

                "Space Filling",

                "Wireframe",

                "Cartoon"
            ]
        )

        show_surface = st.checkbox(
            "Show van der Waals Surface"
        )

        if st.button("🚀 Generate 3D Structure"):

            with st.spinner(
                "Generating 3D molecular structure..."
            ):

                mol3d = generate_3d_molecule(mol)

            if mol3d is None:

                st.error(
                    "Unable to generate the 3D structure."
                )

            else:

                display_3d_structure(

                    mol3d,

                    visualization_style,

                    show_surface
                )


# ============================================================
# COMPLETE MOLECULAR PROPERTIES REPORT
# ============================================================

elif page == "📊 Molecular Properties Report":

    st.title("📊 Complete Molecular Properties Report")

    st.info("""
Enter any valid SMILES notation. The system will automatically
generate a complete molecular properties report.
""")

    smiles = st.text_input(

        "Enter New SMILES",

        value="CCO",

        key="properties_smiles"
    )

    if smiles:

        mol = get_molecule(smiles)

        if mol is None:

            st.error(
                "❌ Invalid SMILES notation. Please check the structure."
            )

        else:

            st.success(
                "✅ Valid Molecular Structure Detected!"
            )

            properties = calculate_properties(mol)

            canonical_smiles = Chem.MolToSmiles(mol)

            # ------------------------------------------------
            # BASIC INFORMATION
            # ------------------------------------------------

            col1, col2 = st.columns([1.2, 1])

            with col1:

                st.subheader("🧬 2D Molecular Structure")

                image = generate_2d_image(mol)

                st.image(image)

            with col2:

                st.subheader("📋 Basic Information")

                st.write(
                    "**Input SMILES:**",
                    smiles
                )

                st.write(
                    "**Canonical SMILES:**",
                    canonical_smiles
                )

                st.write(
                    "**Molecular Formula:**",
                    properties["Molecular Formula"]
                )

                st.write(
                    "**Atoms:**",
                    properties["Number of Atoms"]
                )

                st.write(
                    "**Heavy Atoms:**",
                    properties["Heavy Atoms"]
                )

                st.write(
                    "**Bonds:**",
                    properties["Number of Bonds"]
                )

            st.divider()

            # ------------------------------------------------
            # PROPERTY METRICS
            # ------------------------------------------------

            st.subheader("📊 Important Molecular Properties")

            col1, col2, col3, col4 = st.columns(4)

            col1.metric(
                "Molecular Weight",
                f"{properties['Molecular Weight']} g/mol"
            )

            col2.metric(
                "LogP",
                properties["LogP"]
            )

            col3.metric(
                "TPSA",
                f"{properties['TPSA']} Å²"
            )

            col4.metric(
                "Ring Count",
                properties["Ring Count"]
            )

            col1, col2, col3, col4 = st.columns(4)

            col1.metric(
                "H-Bond Donors",
                properties["Hydrogen Bond Donors"]
            )

            col2.metric(
                "H-Bond Acceptors",
                properties["Hydrogen Bond Acceptors"]
            )

            col3.metric(
                "Rotatable Bonds",
                properties["Rotatable Bonds"]
            )

            col4.metric(
                "Fraction Csp3",
                properties["Fraction Csp3"]
            )

            st.divider()

            # ------------------------------------------------
            # COMPLETE TABLE
            # ------------------------------------------------

            st.subheader("📑 Complete Molecular Property Table")

            report_data = [

                ["Input SMILES", smiles],

                ["Canonical SMILES", canonical_smiles]
            ]

            for key, value in properties.items():

                report_data.append(
                    [key, value]
                )

            report_df = pd.DataFrame(

                report_data,

                columns=[
                    "Property",
                    "Value"
                ]
            )

            st.dataframe(

                report_df,

                use_container_width=True,

                hide_index=True
            )

            # ------------------------------------------------
            # DOWNLOAD CSV
            # ------------------------------------------------

            csv = report_df.to_csv(
                index=False
            ).encode("utf-8")

            st.download_button(

                "⬇️ Download Properties Report (CSV)",

                data=csv,

                file_name="Molecular_Properties_Report.csv",

                mime="text/csv"
            )

            # ------------------------------------------------
            # DOWNLOAD PDF
            # ------------------------------------------------

            if st.button(
                "📄 Generate Molecular PDF Report"
            ):

                with st.spinner(
                    "Generating molecular PDF report..."
                ):

                    pdf = generate_molecule_pdf(

                        smiles,

                        mol,

                        properties
                    )

                st.download_button(

                    "⬇️ Download Molecular Properties PDF",

                    data=pdf,

                    file_name="Molecular_Properties_Report.pdf",

                    mime="application/pdf"
                )

            st.divider()

            # ------------------------------------------------
            # 3D VISUALIZATION
            # ------------------------------------------------

            st.subheader("🌐 3D Molecular Visualization")

            model_style = st.selectbox(

                "Select 3D Model",

                [

                    "Ball and Stick",

                    "Stick",

                    "Space Filling",

                    "Wireframe",

                    "Cartoon"
                ],

                key="properties_3d_style"
            )

            surface = st.checkbox(

                "Show van der Waals Surface",

                key="properties_surface"
            )

            if st.button(
                "🚀 Generate 3D Molecular Model",
                key="properties_generate_3d"
            ):

                with st.spinner(
                    "Generating 3D molecular model..."
                ):

                    mol3d = generate_3d_molecule(mol)

                if mol3d is None:

                    st.error(
                        "Unable to generate the 3D structure."
                    )

                else:

                    display_3d_structure(

                        mol3d,

                        model_style,

                        surface
                    )


# ============================================================
# STRUCTURE PROPERTY ANALYSIS
# ============================================================

elif page == "📈 Structure–Property Analysis":

    st.title("📈 Structure–Property Analysis")

    st.write("""
Compare the properties of multiple molecules.
Enter one molecule per line in the format:

**Molecule Name, SMILES**
""")

    default_data = """Ethanol,CCO
Benzene,c1ccccc1
Phenol,Oc1ccccc1
Acetic Acid,CC(=O)O
Aspirin,CC(=O)Oc1ccccc1C(=O)O"""

    molecule_input = st.text_area(

        "Enter Molecules",

        value=default_data,

        height=200
    )

    rows = []

    for line in molecule_input.splitlines():

        if "," not in line:

            continue

        name, smiles = line.split(",", 1)

        mol = get_molecule(smiles)

        if mol:

            prop = calculate_properties(mol)

            rows.append({

                "Molecule": name.strip(),

                "MW": prop["Molecular Weight"],

                "LogP": prop["LogP"],

                "TPSA": prop["TPSA"],

                "HBD": prop["Hydrogen Bond Donors"],

                "HBA": prop["Hydrogen Bond Acceptors"],

                "Rotatable Bonds": prop["Rotatable Bonds"],

                "Ring Count": prop["Ring Count"]
            })

    if rows:

        df = pd.DataFrame(rows)

        st.dataframe(
            df,
            use_container_width=True,
            hide_index=True
        )

        columns = [

            "MW",

            "LogP",

            "TPSA",

            "HBD",

            "HBA",

            "Rotatable Bonds",

            "Ring Count"
        ]

        col1, col2 = st.columns(2)

        with col1:

            x_axis = st.selectbox(
                "Select X-axis",
                columns
            )

        with col2:

            y_axis = st.selectbox(
                "Select Y-axis",
                columns,
                index=1
            )

        fig, ax = plt.subplots(
            figsize=(8, 5)
        )

        ax.scatter(
            df[x_axis],
            df[y_axis],
            s=100
        )

        for _, row in df.iterrows():

            ax.annotate(

                row["Molecule"],

                (
                    row[x_axis],

                    row[y_axis]
                )
            )

        ax.set_xlabel(x_axis)

        ax.set_ylabel(y_axis)

        ax.set_title(
            f"{x_axis} vs {y_axis}"
        )

        ax.grid(True)

        st.pyplot(fig)

        plt.close(fig)


# ============================================================
# ASSESSMENT
# ============================================================

elif page == "📝 Assessment":

    st.title("📝 Cheminformatics Assessment")

    questions = [

        {
            "question": "What does SMILES represent?",

            "options": [
                "A molecular text representation",
                "A spectroscopy technique",
                "A laboratory instrument"
            ],

            "answer": "A molecular text representation"
        },

        {
            "question": "Which property is commonly associated with molecular lipophilicity?",

            "options": [
                "LogP",
                "TPSA",
                "HBD"
            ],

            "answer": "LogP"
        },

        {
            "question": "What does HBA mean?",

            "options": [
                "Hydrogen Bond Acceptor",
                "Heavy Bond Atom",
                "Hydrogen Bond Analysis"
            ],

            "answer": "Hydrogen Bond Acceptor"
        }
    ]

    student_answers = []

    for i, q in enumerate(questions):

        answer = st.radio(
            q["question"],
            q["options"],
            key=f"assessment_{i}"
        )

        student_answers.append(answer)

    if st.button("Submit Assessment"):

        score = 0

        for i, q in enumerate(questions):

            if student_answers[i] == q["answer"]:
                score += 1

        percentage = (score / len(questions)) * 100

        st.success(
            f"Your Score: {score}/{len(questions)}"
        )

        st.metric(
            "Percentage",
            f"{percentage:.1f}%"
        )

        if percentage >= 80:

            st.balloons()

            st.success("🎉 Excellent performance!")

        elif percentage >= 50:

            st.info("👍 Good! Continue practicing.")

        else:

            st.warning(
                "📚 Please review the Theory section and try again."
            )

# ============================================================
# FINAL PROJECT REPORT
# ============================================================

elif page == "📄 Final Project Report":

    st.title("📄 Final Project Report")

    st.markdown("""
### Cheminformatics Virtual Laboratory

Generate and download the complete final project report
in PDF format.
""")

    if st.button(
        "📄 Generate Final Project Report"
    ):

        with st.spinner(
            "Generating final project report..."
        ):

            pdf = generate_final_report()

        st.success(
            "✅ Final Project Report Generated Successfully!"
        )

        st.download_button(

            label="⬇️ Download Final Project Report PDF",

            data=pdf,

            file_name=(
                "Cheminformatics_Virtual_Laboratory_Final_Report.pdf"
            ),

            mime="application/pdf"
        )


# ============================================================
# FOOTER
# ============================================================

st.divider()

st.caption(
    "🧪 Cheminformatics Virtual Laboratory | "
    "Developed using Python, Streamlit and RDKit"
)
