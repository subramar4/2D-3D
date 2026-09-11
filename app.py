import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from io import BytesIO

from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors

import py3Dmol
import streamlit.components.v1 as components

# PDF GENERATION
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
    PageBreak
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
# HELPER FUNCTIONS
# ============================================================

def get_molecule(smiles):
    """Convert SMILES to RDKit molecule safely."""

    if not smiles or not smiles.strip():
        return None

    try:
        return Chem.MolFromSmiles(smiles.strip())
    except Exception:
        return None


def calculate_properties(mol):
    """Calculate molecular descriptors."""

    return {
        "Molecular Weight": round(Descriptors.MolWt(mol), 2),
        "LogP": round(Crippen.MolLogP(mol), 2),
        "TPSA": round(rdMolDescriptors.CalcTPSA(mol), 2),
        "HBD": Lipinski.NumHDonors(mol),
        "HBA": Lipinski.NumHAcceptors(mol),
        "Rotatable Bonds": Lipinski.NumRotatableBonds(mol),
        "Ring Count": Lipinski.RingCount(mol),
        "Molecular Formula": rdMolDescriptors.CalcMolFormula(mol)
    }


# ============================================================
# PDF REPORT GENERATOR
# ============================================================

def generate_final_report():

    # Create PDF in memory
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
        "CustomTitle",
        parent=styles["Title"],
        fontSize=22,
        leading=28,
        alignment=TA_CENTER,
        textColor=HexColor("#17365D"),
        spaceAfter=15
    )

    subtitle_style = ParagraphStyle(
        "Subtitle",
        parent=styles["Normal"],
        fontSize=13,
        leading=18,
        alignment=TA_CENTER,
        spaceAfter=15
    )

    heading_style = ParagraphStyle(
        "Heading",
        parent=styles["Heading1"],
        fontSize=16,
        leading=20,
        textColor=HexColor("#17365D"),
        spaceBefore=12,
        spaceAfter=8
    )

    subheading_style = ParagraphStyle(
        "SubHeading",
        parent=styles["Heading2"],
        fontSize=12,
        leading=16,
        textColor=HexColor("#2F5597"),
        spaceBefore=8,
        spaceAfter=5
    )

    body_style = ParagraphStyle(
        "Body",
        parent=styles["BodyText"],
        fontSize=10,
        leading=15,
        alignment=TA_JUSTIFY,
        spaceAfter=7
    )

    bullet_style = ParagraphStyle(
        "Bullet",
        parent=body_style,
        leftIndent=18,
        firstLineIndent=-10,
        spaceAfter=4
    )

    story = []

    # ========================================================
    # COVER PAGE
    # ========================================================

    story.append(Spacer(1, 1 * inch))

    story.append(
        Paragraph(
            "FINAL PROJECT REPORT",
            subtitle_style
        )
    )

    story.append(
        Paragraph(
            "Cheminformatics Virtual Laboratory",
            title_style
        )
    )

    story.append(
        Paragraph(
            "A Streamlit-Based Interactive Platform for Molecular "
            "Representation, Visualization, Descriptor Calculation "
            "and Structure–Property Analysis",
            subtitle_style
        )
    )

    story.append(Spacer(1, 0.5 * inch))

    cover_data = [
        ["Project Type", "Interactive Virtual Laboratory"],
        ["Platform", "Streamlit Community Cloud"],
        ["Programming Language", "Python"],
        ["Cheminformatics Toolkit", "RDKit"],
        ["Visualization", "RDKit 2D and py3Dmol 3D"]
    ]

    cover_table = Table(
        cover_data,
        colWidths=[2.1 * inch, 3.5 * inch]
    )

    cover_table.setStyle(
        TableStyle([
            ("BACKGROUND", (0, 0), (0, -1), HexColor("#D9EAF7")),
            ("GRID", (0, 0), (-1, -1), 0.5, colors.grey),
            ("FONTNAME", (0, 0), (0, -1), "Helvetica-Bold"),
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            ("PADDING", (0, 0), (-1, -1), 8)
        ])
    )

    story.append(cover_table)

    story.append(PageBreak())

    # ========================================================
    # ABSTRACT
    # ========================================================

    story.append(
        Paragraph(
            "1. Abstract",
            heading_style
        )
    )

    story.append(
        Paragraph(
            "The Cheminformatics Virtual Laboratory was developed as an "
            "interactive web-based learning environment for introducing "
            "fundamental concepts of computational chemistry and "
            "cheminformatics. The application enables users to study "
            "theoretical concepts, enter molecular structures using SMILES "
            "notation, visualize molecules in two and three dimensions, "
            "calculate important molecular descriptors, compare "
            "structure–property relationships, and complete an assessment. "
            "The complete application was implemented using Python and "
            "Streamlit, with RDKit as the main cheminformatics toolkit.",
            body_style
        )
    )

    # ========================================================
    # OBJECTIVES
    # ========================================================

    story.append(
        Paragraph(
            "2. Objectives",
            heading_style
        )
    )

    objectives = [
        "To develop an interactive virtual laboratory for teaching cheminformatics.",
        "To introduce molecular representation using SMILES notation.",
        "To generate two-dimensional molecular structures.",
        "To visualize molecules in three dimensions.",
        "To calculate important physicochemical molecular descriptors.",
        "To investigate structure–property relationships.",
        "To assess student understanding using an integrated quiz.",
        "To deploy the application using Streamlit Cloud."
    ]

    for item in objectives:

        story.append(
            Paragraph(
                "• " + item,
                bullet_style
            )
        )

    # ========================================================
    # SOFTWARE
    # ========================================================

    story.append(
        Paragraph(
            "3. Software and Technologies Used",
            heading_style
        )
    )

    technology_data = [
        ["Technology", "Purpose"],
        ["Python", "Core programming language."],
        ["Streamlit", "Interactive web application interface."],
        ["RDKit", "Molecular processing and descriptor calculation."],
        ["Pandas", "Data handling and tabular presentation."],
        ["Matplotlib", "Graphical structure–property analysis."],
        ["py3Dmol", "Three-dimensional molecular visualization."],
        ["ReportLab", "PDF report generation."],
        ["Streamlit Cloud", "Application deployment."]
    ]

    technology_table = Table(
        technology_data,
        colWidths=[1.6 * inch, 4.0 * inch],
        repeatRows=1
    )

    technology_table.setStyle(
        TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), HexColor("#17365D")),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
            ("GRID", (0, 0), (-1, -1), 0.5, colors.grey),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("PADDING", (0, 0), (-1, -1), 7)
        ])
    )

    story.append(technology_table)

    # ========================================================
    # MODULES
    # ========================================================

    story.append(
        Paragraph(
            "4. Functional Modules",
            heading_style
        )
    )

    modules = [

        (
            "4.1 Home Module",
            "Introduces the virtual laboratory, learning objectives "
            "and example molecules."
        ),

        (
            "4.2 Theory Module",
            "Provides theoretical information about cheminformatics, "
            "SMILES notation and molecular descriptors."
        ),

        (
            "4.3 Molecular Visualization",
            "Allows users to enter SMILES notation and generate "
            "two-dimensional molecular structures."
        ),

        (
            "4.4 Three-Dimensional Visualization",
            "Generates molecular geometry and displays interactive "
            "three-dimensional structures."
        ),

        (
            "4.5 Molecular Descriptor Calculator",
            "Calculates Molecular Weight, LogP, TPSA, HBD, HBA, "
            "Rotatable Bonds and Ring Count."
        ),

        (
            "4.6 Structure–Property Analysis",
            "Compares molecular properties and generates graphical "
            "relationships between selected descriptors."
        ),

        (
            "4.7 Assessment Module",
            "Evaluates student understanding using multiple-choice "
            "questions and automatic scoring."
        )
    ]

    for heading, description in modules:

        story.append(
            Paragraph(
                heading,
                subheading_style
            )
        )

        story.append(
            Paragraph(
                description,
                body_style
            )
        )

    story.append(PageBreak())

    # ========================================================
    # MOLECULAR DESCRIPTORS
    # ========================================================

    story.append(
        Paragraph(
            "5. Molecular Descriptors",
            heading_style
        )
    )

    descriptor_data = [
        ["Descriptor", "Description"],
        ["Molecular Weight", "Total molecular mass."],
        ["LogP", "Measure associated with molecular lipophilicity."],
        ["TPSA", "Topological polar surface area."],
        ["HBD", "Number of hydrogen bond donors."],
        ["HBA", "Number of hydrogen bond acceptors."],
        ["Rotatable Bonds", "Indicator of molecular flexibility."],
        ["Ring Count", "Number of molecular rings."]
    ]

    descriptor_table = Table(
        descriptor_data,
        colWidths=[1.7 * inch, 3.9 * inch],
        repeatRows=1
    )

    descriptor_table.setStyle(
        TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), HexColor("#2F5597")),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
            ("GRID", (0, 0), (-1, -1), 0.5, colors.grey),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("PADDING", (0, 0), (-1, -1), 7)
        ])
    )

    story.append(descriptor_table)

    # ========================================================
    # DEPLOYMENT
    # ========================================================

    story.append(
        Paragraph(
            "6. Deployment and Dependency Configuration",
            heading_style
        )
    )

    story.append(
        Paragraph(
            "The application was deployed using Streamlit Community Cloud. "
            "During deployment, dependency-related errors were encountered "
            "and resolved through proper configuration of Python packages "
            "and Linux system libraries.",
            body_style
        )
    )

    story.append(
        Paragraph(
            "6.1 Python Dependencies",
            subheading_style
        )
    )

    requirements_data = [
        ["Package"],
        ["streamlit"],
        ["pandas"],
        ["matplotlib"],
        ["py3Dmol"],
        ["rdkit"],
        ["reportlab"]
    ]

    requirements_table = Table(
        requirements_data,
        colWidths=[3 * inch]
    )

    requirements_table.setStyle(
        TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), HexColor("#17365D")),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
            ("GRID", (0, 0), (-1, -1), 0.5, colors.grey),
            ("ALIGN", (0, 0), (-1, -1), "CENTER"),
            ("PADDING", (0, 0), (-1, -1), 6)
        ])
    )

    story.append(requirements_table)

    story.append(
        Paragraph(
            "6.2 Linux System Dependencies",
            subheading_style
        )
    )

    story.append(
        Paragraph(
            "The RDKit drawing functionality required additional Linux "
            "libraries. The missing libXrender.so.1 error was resolved "
            "using packages.txt.",
            body_style
        )
    )

    packages_data = [
        ["System Package"],
        ["libxrender1"],
        ["libxext6"],
        ["libsm6"]
    ]

    packages_table = Table(
        packages_data,
        colWidths=[3 * inch]
    )

    packages_table.setStyle(
        TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), HexColor("#17365D")),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
            ("GRID", (0, 0), (-1, -1), 0.5, colors.grey),
            ("ALIGN", (0, 0), (-1, -1), "CENTER"),
            ("PADDING", (0, 0), (-1, -1), 6)
        ])
    )

    story.append(packages_table)

    # ========================================================
    # PROBLEMS AND SOLUTIONS
    # ========================================================

    story.append(
        Paragraph(
            "7. Problems Encountered and Solutions",
            heading_style
        )
    )

    problem_data = [
        ["Problem", "Solution"],

        [
            "Missing Python packages",
            "Added the required packages to requirements.txt."
        ],

        [
            "libXrender.so.1 missing",
            "Added libxrender1, libxext6 and libsm6 to packages.txt."
        ],

        [
            "PDF file not found",
            "Generated the PDF directly in memory using BytesIO."
        ],

        [
            "Single-file Streamlit application",
            "Used sidebar navigation and one page configuration."
        ]
    ]

    problem_table = Table(
        problem_data,
        colWidths=[2.2 * inch, 3.4 * inch],
        repeatRows=1
    )

    problem_table.setStyle(
        TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), HexColor("#2F5597")),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
            ("GRID", (0, 0), (-1, -1), 0.5, colors.grey),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("PADDING", (0, 0), (-1, -1), 7)
        ])
    )

    story.append(problem_table)

    # ========================================================
    # FINAL OUTCOME
    # ========================================================

    story.append(
        Paragraph(
            "8. Final Outcome",
            heading_style
        )
    )

    final_points = [

        "Interactive single-file Streamlit application.",

        "SMILES input and molecular validation.",

        "Two-dimensional molecular visualization.",

        "Three-dimensional molecular visualization.",

        "Automatic molecular descriptor calculation.",

        "Structure–property comparison and graphs.",

        "Integrated student assessment.",

        "PDF final report generation and download.",

        "Successful deployment configuration."
    ]

    for item in final_points:

        story.append(
            Paragraph(
                "• " + item,
                bullet_style
            )
        )

    # ========================================================
    # FUTURE SCOPE
    # ========================================================

    story.append(
        Paragraph(
            "9. Future Scope",
            heading_style
        )
    )

    future_scope = [

        "Addition of MOL and SDF file upload.",

        "More advanced molecular descriptors.",

        "Machine learning-based property prediction.",

        "Additional student assessments.",

        "Student progress tracking.",

        "More advanced molecular visualization tools."
    ]

    for item in future_scope:

        story.append(
            Paragraph(
                "• " + item,
                bullet_style
            )
        )

    # ========================================================
    # CONCLUSION
    # ========================================================

    story.append(
        Paragraph(
            "10. Conclusion",
            heading_style
        )
    )

    story.append(
        Paragraph(
            "The Cheminformatics Virtual Laboratory successfully integrates "
            "Python-based scientific computing with interactive web technology "
            "for chemical education. The combination of Streamlit and RDKit "
            "provides an accessible environment for molecular representation, "
            "visualization, descriptor calculation and structure–property "
            "analysis. The final application also includes PDF report generation, "
            "allowing users to download the project report directly from the "
            "web application.",
            body_style
        )
    )

    # ========================================================
    # FOOTER
    # ========================================================

    def add_page_number(canvas, doc):

        canvas.saveState()

        canvas.setFont("Helvetica", 8)

        canvas.setFillColor(colors.grey)

        canvas.drawCentredString(
            A4[0] / 2,
            25,
            f"Cheminformatics Virtual Laboratory | Page {doc.page}"
        )

        canvas.restoreState()

    # Build PDF
    doc.build(
        story,
        onFirstPage=add_page_number,
        onLaterPages=add_page_number
    )

    # IMPORTANT: Get PDF bytes from memory
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
        "🧬 SMILES & Molecular Visualization",
        "📊 Molecular Descriptors",
        "📈 Structure–Property Analysis",
        "📝 Assessment",
        "📄 Final Report"
    ]
)


st.sidebar.divider()

st.sidebar.info("""
### Learning Objectives

• SMILES representation

• 2D and 3D visualization

• Molecular descriptors

• Structure–property relationships
""")


# ============================================================
# HOME
# ============================================================

if page == "🏠 Home":

    st.title("🧪 Cheminformatics Virtual Laboratory")

    st.subheader(
        "Computational Representation, Visualization and Analysis "
        "of Chemical Structures"
    )

    st.markdown("""
## Welcome

This virtual laboratory introduces fundamental concepts of
**Cheminformatics**.

Students can:

- Learn molecular representation
- Enter SMILES notation
- Visualize molecules
- Calculate descriptors
- Analyze structure–property relationships
- Complete an assessment
""")

    st.markdown("## 🧪 Example Molecules")

    examples = pd.DataFrame({

        "Molecule": [
            "Ethanol",
            "Benzene",
            "Acetic Acid",
            "Caffeine"
        ],

        "SMILES": [
            "CCO",
            "c1ccccc1",
            "CC(=O)O",
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
        ]
    })

    st.dataframe(
        examples,
        use_container_width=True,
        hide_index=True
    )


# ============================================================
# THEORY
# ============================================================

elif page == "📖 Theory":

    st.title("📖 Introduction to Cheminformatics")

    st.markdown("""
## What is Cheminformatics?

Cheminformatics combines:

- ⚗️ Chemistry
- 💻 Computer Science
- 📊 Data Analysis

It uses computational methods to represent, store and analyze
chemical information.
""")

    st.markdown("""
## SMILES

SMILES means:

**Simplified Molecular Input Line Entry System**

It represents molecular structures using text.

Examples:

- Ethanol → `CCO`
- Benzene → `c1ccccc1`
- Acetic Acid → `CC(=O)O`
""")

    st.markdown("""
## Molecular Descriptors

Important descriptors include:

- Molecular Weight
- LogP
- TPSA
- Hydrogen Bond Donors
- Hydrogen Bond Acceptors
""")


# ============================================================
# VISUALIZATION
# ============================================================

elif page == "🧬 SMILES & Molecular Visualization":

    st.title("🧬 Molecular Visualization")

    examples_dict = {

        "Ethanol": "CCO",

        "Benzene": "c1ccccc1",

        "Acetic Acid": "CC(=O)O",

        "Caffeine":
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    }

    if "visualization_smiles" not in st.session_state:

        st.session_state.visualization_smiles = "CCO"

    col1, col2 = st.columns([2, 1])

    with col1:

        st.text_input(
            "Enter SMILES",
            key="visualization_smiles"
        )

    with col2:

        selected_example = st.selectbox(
            "Choose Example",
            list(examples_dict.keys())
        )

        if st.button("Load Example"):

            st.session_state.visualization_smiles = (
                examples_dict[selected_example]
            )

            st.rerun()

    smiles = st.session_state.visualization_smiles

    mol = get_molecule(smiles)

    if mol is None:

        st.error("❌ Invalid SMILES")

    else:

        st.success("✅ Valid molecular structure")

        col1, col2 = st.columns(2)

        with col1:

            st.subheader("2D Structure")

            image = Draw.MolToImage(
                mol,
                size=(500, 400)
            )

            st.image(image)

        with col2:

            properties = calculate_properties(mol)

            st.subheader("Molecular Information")

            st.write(
                "**Canonical SMILES:**",
                Chem.MolToSmiles(mol)
            )

            st.write(
                "**Molecular Formula:**",
                properties["Molecular Formula"]
            )

            st.write(
                "**Number of Atoms:**",
                mol.GetNumAtoms()
            )

        st.divider()

        if st.button("🌐 Generate 3D Structure"):

            try:

                mol3d = Chem.AddHs(mol)

                status = AllChem.EmbedMolecule(
                    mol3d,
                    randomSeed=42
                )

                if status == 0:

                    try:
                        AllChem.MMFFOptimizeMolecule(mol3d)
                    except Exception:
                        pass

                    mol_block = Chem.MolToMolBlock(mol3d)

                    view = py3Dmol.view(
                        width=900,
                        height=500
                    )

                    view.addModel(
                        mol_block,
                        "mol"
                    )

                    view.setStyle({
                        "stick": {}
                    })

                    view.setBackgroundColor("white")

                    view.zoomTo()

                    components.html(
                        view._make_html(),
                        height=520
                    )

                else:

                    st.warning(
                        "Unable to generate 3D structure."
                    )

            except Exception as e:

                st.error(f"3D Error: {e}")


# ============================================================
# MOLECULAR DESCRIPTORS
# ============================================================

elif page == "📊 Molecular Descriptors":

    st.title("📊 Molecular Descriptor Calculator")

    smiles = st.text_input(
        "Enter SMILES",
        value="CCO",
        key="descriptor_smiles"
    )

    mol = get_molecule(smiles)

    if mol is None:

        st.error("❌ Invalid SMILES")

    else:

        properties = calculate_properties(mol)

        col1, col2, col3, col4, col5 = st.columns(5)

        col1.metric(
            "Molecular Weight",
            properties["Molecular Weight"]
        )

        col2.metric(
            "LogP",
            properties["LogP"]
        )

        col3.metric(
            "TPSA",
            properties["TPSA"]
        )

        col4.metric(
            "HBD",
            properties["HBD"]
        )

        col5.metric(
            "HBA",
            properties["HBA"]
        )

        descriptor_df = pd.DataFrame({

            "Descriptor": list(properties.keys()),

            "Value": list(properties.values())
        })

        st.dataframe(
            descriptor_df,
            use_container_width=True,
            hide_index=True
        )


# ============================================================
# STRUCTURE-PROPERTY ANALYSIS
# ============================================================

elif page == "📈 Structure–Property Analysis":

    st.title("📈 Structure–Property Analysis")

    default_data = """Ethanol,CCO
Benzene,c1ccccc1
Phenol,Oc1ccccc1
Acetic Acid,CC(=O)O
Caffeine,CN1C=NC2=C1C(=O)N(C(=O)N2C)C"""

    molecule_input = st.text_area(
        "Enter: Name,SMILES",
        value=default_data,
        height=200
    )

    rows = []

    for line in molecule_input.splitlines():

        if "," not in line:
            continue

        name, smiles = line.split(",", 1)

        mol = get_molecule(smiles)

        if mol is not None:

            prop = calculate_properties(mol)

            rows.append({

                "Molecule": name.strip(),

                "MW": prop["Molecular Weight"],

                "LogP": prop["LogP"],

                "TPSA": prop["TPSA"],

                "HBD": prop["HBD"],

                "HBA": prop["HBA"]
            })

    if rows:

        df = pd.DataFrame(rows)

        st.dataframe(
            df,
            use_container_width=True,
            hide_index=True
        )

        properties_list = [
            "MW",
            "LogP",
            "TPSA",
            "HBD",
            "HBA"
        ]

        col1, col2 = st.columns(2)

        with col1:

            x_axis = st.selectbox(
                "X-axis",
                properties_list,
                index=0
            )

        with col2:

            y_axis = st.selectbox(
                "Y-axis",
                properties_list,
                index=1
            )

        fig, ax = plt.subplots()

        ax.scatter(
            df[x_axis],
            df[y_axis],
            s=100
        )

        for _, row in df.iterrows():

            ax.annotate(
                row["Molecule"],
                (row[x_axis], row[y_axis])
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

    st.title("📝 Assessment")

    questions = [

        (
            "What does SMILES represent?",

            [
                "A molecular text representation",
                "A spectroscopy technique",
                "A laboratory instrument"
            ],

            "A molecular text representation"
        ),

        (
            "Which descriptor represents lipophilicity?",

            [
                "TPSA",
                "LogP",
                "HBD"
            ],

            "LogP"
        ),

        (
            "HBD means?",

            [
                "Hydrogen Bond Donor",
                "High Bond Density",
                "Heavy Bond Data"
            ],

            "Hydrogen Bond Donor"
        )
    ]

    answers = []

    for i, question in enumerate(questions):

        answer = st.radio(
            question[0],
            question[1],
            key=f"question_{i}"
        )

        answers.append(answer)

    if st.button("Submit Assessment"):

        score = 0

        for i, question in enumerate(questions):

            if answers[i] == question[2]:

                score += 1

        percentage = (
            score / len(questions)
        ) * 100

        st.success(
            f"Score: {score}/{len(questions)}"
        )

        st.metric(
            "Percentage",
            f"{percentage:.1f}%"
        )


# ============================================================
# FINAL REPORT
# ============================================================

elif page == "📄 Final Report":

    st.title("📄 Final Project Report")

    st.markdown("""
### Cheminformatics Virtual Laboratory

You can generate and download the complete final project report
in PDF format.
""")

    if st.button("📄 Generate Final Report PDF"):

        with st.spinner("Generating PDF report..."):

            pdf = generate_final_report()

        st.success("✅ PDF Report Generated Successfully!")

        st.download_button(

            label="⬇️ Download Final Report PDF",

            data=pdf,

            file_name=(
                "Cheminformatics_Virtual_Laboratory_Final_Report.pdf"
            ),

            mime="application/pdf"
        )
