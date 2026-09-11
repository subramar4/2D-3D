import streamlit as st
import pandas as pd
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
from reportlab.lib.utils import ImageReader
import tempfile
import os


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
# CALCULATE SIGMA AND PI BONDS
# ============================================================

def calculate_bond_information(mol):
    """
    Calculate:
    - Heavy-atom bonds
    - Sigma (σ) bonds, including bonds to implicit hydrogens
    - Pi (π) bonds
    - Total bond order count (σ + π)
    """

    # Heavy-atom graph bonds (RDKit default molecule representation)
    heavy_atom_bonds = mol.GetNumBonds()

    # Count π bonds using a Kekulé representation for aromatic systems
    kekule_mol = Chem.Mol(mol)

    try:
        Chem.Kekulize(kekule_mol, clearAromaticFlags=True)
    except Exception:
        # If Kekulization is not possible, continue with the available form
        pass

    pi_bonds = 0

    for bond in kekule_mol.GetBonds():
        bond_type = bond.GetBondType()

        if bond_type == Chem.BondType.DOUBLE:
            pi_bonds += 1
        elif bond_type == Chem.BondType.TRIPLE:
            pi_bonds += 2
        elif bond_type == Chem.BondType.QUADRUPLE:
            pi_bonds += 3
        elif bond.GetIsAromatic():
            # Fallback for aromatic bonds if Kekulization did not clear them
            pi_bonds += 0.5

    # Add implicit hydrogens. Every atom-to-atom connection contains one σ bond.
    mol_with_h = Chem.AddHs(Chem.Mol(mol))
    sigma_bonds = mol_with_h.GetNumBonds()

    total_bonds = sigma_bonds + pi_bonds

    return {
        "Heavy Atom Bonds": heavy_atom_bonds,
        "Sigma Bonds": sigma_bonds,
        "Pi Bonds": pi_bonds,
        "Total Bonds": total_bonds,
        "Total Atoms Including H": mol_with_h.GetNumAtoms(),
    }


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
# REPORT IMAGE HELPERS
# ============================================================

def save_2d_structure_png(mol):
    """Create a temporary PNG for inclusion in PDF reports."""
    img = Draw.MolToImage(mol, size=(700, 500))
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".png")
    tmp.close()
    img.save(tmp.name)
    return tmp.name


def save_3d_structure_png(mol):
    """Create a simple 3D coordinate rendering from an RDKit conformer."""
    mol3d = generate_3d_molecule(mol)
    if mol3d is None:
        return None

    conf = mol3d.GetConformer()
    xs, ys, zs, labels = [], [], [], []
    for atom in mol3d.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        xs.append(pos.x); ys.append(pos.y); zs.append(pos.z)
        labels.append(atom.GetSymbol())

    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111, projection="3d")

    for bond in mol3d.GetBonds():
        a = conf.GetAtomPosition(bond.GetBeginAtomIdx())
        b = conf.GetAtomPosition(bond.GetEndAtomIdx())
        ax.plot([a.x, b.x], [a.y, b.y], [a.z, b.z])

    ax.scatter(xs, ys, zs, s=80)
    for x, y, z, label in zip(xs, ys, zs, labels):
        ax.text(x, y, z, label, fontsize=8)

    ax.set_xlabel("X"); ax.set_ylabel("Y"); ax.set_zlabel("Z")
    ax.set_title("3D Molecular Structure")
    ax.view_init(elev=25, azim=45)

    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".png")
    tmp.close()
    fig.savefig(tmp.name, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return tmp.name


# ============================================================
# GENERATE GENERAL FINAL PROJECT REPORT
# ============================================================

def generate_final_report(student_name, registration_number, input_smiles, mol, properties):
    """Generate a personalized final report for the student's entered molecule."""
    buffer = BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=A4, leftMargin=45, rightMargin=45,
                            topMargin=45, bottomMargin=45)
    styles = getSampleStyleSheet()
    title_style = ParagraphStyle("ProjectTitle", parent=styles["Title"], fontSize=20,
                                 leading=25, alignment=TA_CENTER, textColor=HexColor("#17365D"))
    heading_style = ParagraphStyle("ProjectHeading", parent=styles["Heading1"], fontSize=14,
                                   leading=18, textColor=HexColor("#17365D"))
    body_style = ParagraphStyle("ProjectBody", parent=styles["BodyText"], fontSize=10,
                                leading=14, alignment=TA_JUSTIFY)

    story = []
    story += [Spacer(1, 0.5*inch), Paragraph("CHEMINFORMATICS VIRTUAL LABORATORY", title_style),
              Spacer(1, 0.2*inch), Paragraph("FINAL STUDENT MOLECULAR ANALYSIS REPORT", heading_style),
              Spacer(1, 0.25*inch)]

    student_data = [["Student Name", student_name or "Not provided"],
                    ["Registration Number", registration_number or "Not provided"],
                    ["Input SMILES", input_smiles],
                    ["Canonical SMILES", Chem.MolToSmiles(mol)],
                    ["Molecular Formula", properties["Molecular Formula"]]]
    t = Table(student_data, colWidths=[2.1*inch, 4.0*inch])
    t.setStyle(TableStyle([("GRID", (0,0),(-1,-1),0.5,colors.grey),
                           ("BACKGROUND",(0,0),(0,-1),HexColor("#EAF0F8")),
                           ("VALIGN",(0,0),(-1,-1),"TOP"),
                           ("PADDING",(0,0),(-1,-1),6)]))
    story += [t, Spacer(1,0.25*inch)]

    story.append(Paragraph("1. Molecular Structures", heading_style))
    img2d = save_2d_structure_png(mol)
    img3d = save_3d_structure_png(mol)
    try:
        if img2d and img3d:
            structure_table = Table([[Image(img2d, width=3.0*inch, height=2.2*inch),
                                      Image(img3d, width=3.0*inch, height=2.2*inch)]],
                                    colWidths=[3.1*inch,3.1*inch])
            structure_table.setStyle(TableStyle([("VALIGN",(0,0),(-1,-1),"MIDDLE")]))
            story += [structure_table, Paragraph("Left: 2D molecular structure &nbsp;&nbsp;&nbsp; Right: generated 3D molecular structure", body_style)]
        elif img2d:
            story.append(Image(img2d, width=5.5*inch, height=3.9*inch))
    finally:
        pass
    story += [Spacer(1,0.2*inch)]

    story.append(Paragraph("2. Complete Molecular Properties", heading_style))
    prop_data = [["Property", "Value"]]
    for k,v in properties.items(): prop_data.append([k, str(v)])
    pt = Table(prop_data, colWidths=[3.5*inch,2.7*inch], repeatRows=1)
    pt.setStyle(TableStyle([("BACKGROUND",(0,0),(-1,0),HexColor("#2F5597")),
                            ("TEXTCOLOR",(0,0),(-1,0),colors.white),
                            ("GRID",(0,0),(-1,-1),0.5,colors.grey),
                            ("VALIGN",(0,0),(-1,-1),"TOP"),
                            ("PADDING",(0,0),(-1,-1),5)]))
    story += [pt, Spacer(1,0.2*inch)]

    story.append(Paragraph("3. Interpretation", heading_style))
    story.append(Paragraph(
        "This personalized report was generated from the student's submitted SMILES notation. "
        "The application converted the molecular representation into a validated chemical structure, "
        "calculated physicochemical descriptors and bond information, and generated both 2D and 3D molecular representations using RDKit.",
        body_style))

    story += [Spacer(1,0.2*inch), Paragraph("Developed by Dr. R. Subramanian, Assistant Professor, Division of Chemistry, School of Sciences, Faculty of Engineering and Technology, SRM Institute of Science and Technology, Tiruchirappalli.", body_style)]
    story.append(Paragraph("For queries: rsmani84@gmail.com", body_style))

    doc.build(story)
    pdf_data = buffer.getvalue(); buffer.close()
    for path in [img2d, img3d]:
        if path and os.path.exists(path):
            try: os.remove(path)
            except Exception: pass
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

st.sidebar.markdown("### 👩‍🎓 Student Information")
student_name = st.sidebar.text_input("Student Name", key="student_name")
registration_number = st.sidebar.text_input("Registration Number", key="registration_number")

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

    st.divider()

    st.markdown("### 👨‍🏫 Developed by")
    st.markdown("**Dr. R. Subramanian**")
    st.markdown("Assistant Professor  ")
    st.markdown("Division of Chemistry, School of Sciences  ")
    st.markdown("Faculty of Engineering and Technology  ")
    st.markdown("SRM Institute of Science and Technology, Tiruchirappalli")

    st.info("📧 **For queries, contact:** rsmani84@gmail.com")


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

    # Initialize widget state BEFORE the text_input is created.
    # This avoids StreamlitWidgetAlreadyInstantiatedError.
    if "visualization_smiles" not in st.session_state:
        st.session_state.visualization_smiles = "CCO"

    if "selected_example" not in st.session_state:
        st.session_state.selected_example = "Ethanol"

    def load_selected_example():
        st.session_state.visualization_smiles = (
            examples[st.session_state.selected_example]
        )

    col1, col2 = st.columns([2, 1])

    with col1:
        smiles = st.text_input(
            "Enter SMILES",
            key="visualization_smiles"
        )

    with col2:
        st.selectbox(
            "Choose Example Molecule",
            list(examples.keys()),
            key="selected_example"
        )

        st.button(
            "Load Example",
            on_click=load_selected_example,
            key="load_visualization_example"
        )

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

            bond_info = calculate_bond_information(mol)

            st.write(
                "**Heavy Atoms:**",
                mol.GetNumHeavyAtoms()
            )

            st.write(
                "**Total Atoms (including H):**",
                bond_info["Total Atoms Including H"]
            )

            st.write(
                "**Heavy-Atom Bonds:**",
                bond_info["Heavy Atom Bonds"]
            )

            st.write(
                "**Sigma (σ) Bonds:**",
                bond_info["Sigma Bonds"]
            )

            st.write(
                "**Pi (π) Bonds:**",
                bond_info["Pi Bonds"]
            )

            st.write(
                "**Total Bonds (σ + π):**",
                bond_info["Total Bonds"]
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

    st.markdown("### 🧪 Enter a Molecule for Analysis")
    st.caption("Enter a new SMILES below. The complete report will be generated automatically for that molecule.")

    smiles = st.text_input(
        "Enter New SMILES",
        placeholder="Example: CCO or Oc1ccccc1",
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

            # Save the latest student-entered molecule for the Final Project Report
            st.session_state.latest_smiles = smiles
            st.session_state.latest_properties = properties

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
                    "**Heavy-Atom Bonds:**",
                    properties["Heavy Atom Bonds"]
                )

                st.write(
                    "**Sigma (σ) Bonds:**",
                    properties["Sigma (σ) Bonds"]
                )

                st.write(
                    "**Pi (π) Bonds:**",
                    properties["Pi (π) Bonds"]
                )

                st.write(
                    "**Total Bonds (σ + π):**",
                    properties["Total Bonds (σ + π)"]
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

        st.info("📋 The comparison graph has been removed. The table above provides the molecular property comparison directly.")

        csv_comparison = df.to_csv(index=False).encode("utf-8")
        st.download_button(
            "⬇️ Download Comparison Table (CSV)",
            data=csv_comparison,
            file_name="Structure_Property_Comparison.csv",
            mime="text/csv"
        )

# ============================================================
# ASSESSMENT
# ============================================================

elif page == "📝 Assessment":

    st.title("📝 Cheminformatics Assessment")
    st.info("Select your answers and submit the assessment. Only your score and percentage will be displayed; correct answers will not be revealed.")

    questions = [
        {"question":"What is the full form of SMILES?", "options":["Simplified Molecular Input Line Entry System","Standard Molecular Information Language Encoding System","Simple Molecular Identification and Labeling System"], "answer":"Simplified Molecular Input Line Entry System"},
        {"question":"What does a SMILES string represent?", "options":["A molecular structure using text","A spectroscopy technique","A laboratory instrument"], "answer":"A molecular structure using text"},
        {"question":"Which SMILES represents ethanol?", "options":["CCO","c1ccccc1","CC(=O)O"], "answer":"CCO"},
        {"question":"Which SMILES represents benzene?", "options":["c1ccccc1","CCO","O=C=O"], "answer":"c1ccccc1"},
        {"question":"In SMILES, parentheses are mainly used to indicate:", "options":["Branches","Molecular weight","Atom colour"], "answer":"Branches"},
        {"question":"Which property is commonly associated with molecular lipophilicity?", "options":["LogP","TPSA","HBD"], "answer":"LogP"},
        {"question":"What does TPSA describe?", "options":["Topological Polar Surface Area","Total Pi Surface Area","Thermal Property Surface Analysis"], "answer":"Topological Polar Surface Area"},
        {"question":"What does HBA mean?", "options":["Hydrogen Bond Acceptor","Heavy Bond Atom","Hydrogen Bond Analysis"], "answer":"Hydrogen Bond Acceptor"},
        {"question":"What does HBD mean?", "options":["Hydrogen Bond Donor","Heavy Bond Descriptor","Hydrogen Bond Distance"], "answer":"Hydrogen Bond Donor"},
        {"question":"A double bond contains how many pi (π) bonds?", "options":["1","2","0"], "answer":"1"},
        {"question":"A triple bond contains how many pi (π) bonds?", "options":["2","1","3"], "answer":"2"},
        {"question":"Which model is useful for visualizing approximate molecular volume?", "options":["Space Filling","Wireframe only","Text SMILES"], "answer":"Space Filling"}
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

    st.title("📄 Final Student Molecular Analysis Report")
    st.write("This report uses the latest valid SMILES entered in the Molecular Properties Report page.")

    name = st.text_input("Student Name for Report", value=st.session_state.get("student_name", ""), key="final_report_name")
    reg = st.text_input("Registration Number for Report", value=st.session_state.get("registration_number", ""), key="final_report_reg")

    default_smiles = st.session_state.get("latest_smiles", "CCO")
    report_smiles = st.text_input("SMILES for Final Report", value=default_smiles, key="final_report_smiles")
    report_mol = get_molecule(report_smiles)

    if report_mol is None:
        st.error("❌ Please enter a valid SMILES notation for the final report.")
    else:
        report_properties = calculate_properties(report_mol)
        st.success("✅ Valid molecule ready for the final report.")
        st.write("**Molecular Formula:**", report_properties["Molecular Formula"])
        st.write("**Molecular Weight:**", report_properties["Molecular Weight"], "g/mol")
        st.write("**Sigma Bonds:**", report_properties["Sigma (σ) Bonds"])
        st.write("**Pi Bonds:**", report_properties["Pi (π) Bonds"])

        if st.button("📄 Generate Personalized Final Project Report"):
            with st.spinner("Generating your personalized report with molecular structures and properties..."):
                pdf = generate_final_report(name, reg, report_smiles, report_mol, report_properties)
            st.success("✅ Final Project Report Generated Successfully!")
            st.download_button(
                label="⬇️ Download Personalized Final Project Report PDF",
                data=pdf,
                file_name="Cheminformatics_Final_Molecular_Report.pdf",
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
