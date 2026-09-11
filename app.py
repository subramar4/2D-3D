import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt

from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors

import py3Dmol
import streamlit.components.v1 as components


# ============================================================
# PAGE CONFIGURATION
# IMPORTANT: Use ONLY ONCE
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
    return Chem.MolFromSmiles(smiles.strip())


def calculate_properties(mol):
    """Calculate important molecular descriptors."""

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
# SIDEBAR NAVIGATION
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
        "📝 Assessment"
    ],
    key="main_navigation"
)

st.sidebar.divider()

st.sidebar.info("""
### Learning Objectives

• SMILES representation

• 2D and 3D visualization

• Molecular descriptors

• Structure–property relationships
""")

st.sidebar.divider()

st.sidebar.caption("🧪 Cheminformatics Virtual Laboratory")


# ============================================================
# PAGE 1: HOME
# ============================================================

if page == "🏠 Home":

    st.title("🧪 Cheminformatics Virtual Laboratory")

    st.subheader(
        "Computational Representation, Visualization and Analysis "
        "of Chemical Structures"
    )

    st.divider()

    st.markdown("""
## Welcome to the Virtual Laboratory

This experiment introduces students to the fundamental concepts of
**Cheminformatics** using computational tools.

Students will learn how to represent chemical structures using
**SMILES notation**, visualize molecules in **2D and 3D**, calculate
important molecular descriptors, and investigate
**structure–property relationships**.
""")

    st.markdown("## 🎯 Learning Objectives")

    col1, col2 = st.columns(2)

    with col1:

        st.success("""
### 1. Molecular Representation

Convert molecular representations using **SMILES notation**.
""")

        st.success("""
### 2. Molecular Visualization

Visualize chemical structures in **2D and 3D**.
""")

    with col2:

        st.success("""
### 3. Molecular Descriptor Calculation

Calculate Molecular Weight, LogP, TPSA, HBD and HBA.
""")

        st.success("""
### 4. Structure–Property Relationships

Understand how molecular structure influences chemical properties.
""")

    st.divider()

    st.markdown("## 🔬 Virtual Laboratory Workflow")

    col1, col2, col3, col4 = st.columns(4)

    with col1:
        st.info("### Step 1\n\n📖 Learn Theory")

    with col2:
        st.info("### Step 2\n\n🧬 Enter SMILES")

    with col3:
        st.info("### Step 3\n\n📊 Calculate Properties")

    with col4:
        st.info("### Step 4\n\n📈 Analyze Results")

    st.divider()

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
# PAGE 2: THEORY
# ============================================================

elif page == "📖 Theory":

    st.title("📖 Theory: Introduction to Cheminformatics")

    st.markdown("""
## What is Cheminformatics?

**Cheminformatics** is a scientific field that combines:

- Chemistry ⚗️
- Computer Science 💻
- Data Analysis 📊

It uses computational methods to store, represent, visualize,
analyze and predict chemical information.
""")

    st.divider()

    st.markdown("""
## 🧬 Molecular Representation

Molecules can be represented in several ways.

### Molecular Formula

Shows the number and type of atoms.

Example:

**Ethanol → C₂H₆O**

### Structural Formula

Shows how atoms are connected.

Example:

**CH₃–CH₂–OH**

### SMILES Notation

SMILES stands for:

**Simplified Molecular Input Line Entry System**

It represents a molecular structure as text.
""")

    examples = pd.DataFrame({
        "Molecule": [
            "Water",
            "Ethanol",
            "Benzene",
            "Acetic Acid"
        ],
        "SMILES": [
            "O",
            "CCO",
            "c1ccccc1",
            "CC(=O)O"
        ]
    })

    st.dataframe(
        examples,
        use_container_width=True,
        hide_index=True
    )

    st.divider()

    st.markdown("""
## 📊 Molecular Descriptors

Molecular descriptors are numerical values that describe
important molecular properties.
""")

    descriptor_table = pd.DataFrame({
        "Descriptor": [
            "Molecular Weight",
            "LogP",
            "TPSA",
            "HBD",
            "HBA"
        ],
        "Meaning": [
            "Molecular mass",
            "Lipophilicity",
            "Molecular polarity",
            "Hydrogen bond donors",
            "Hydrogen bond acceptors"
        ]
    })

    st.table(descriptor_table)

    st.markdown("""
## 🔬 Structure–Property Relationship

Changes in molecular structure can influence:

- Molecular weight
- Polarity
- Solubility
- Lipophilicity
- Hydrogen bonding ability

Functional groups and molecular size are important factors
controlling these properties.
""")


# ============================================================
# PAGE 3: MOLECULAR VISUALIZATION
# ============================================================

elif page == "🧬 SMILES & Molecular Visualization":

    st.title("🧬 SMILES and Molecular Visualization")

    st.markdown("""
Enter a valid **SMILES string** and convert it into a
chemical structure.
""")

    examples_dict = {
        "Ethanol": "CCO",
        "Benzene": "c1ccccc1",
        "Acetic Acid": "CC(=O)O",
        "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    }

    # Initialize session state
    if "visualization_smiles" not in st.session_state:
        st.session_state.visualization_smiles = "CCO"

    col1, col2 = st.columns([2, 1])

    with col1:

        smiles = st.text_input(
            "Enter SMILES",
            key="visualization_smiles"
        )

    with col2:

        example = st.selectbox(
            "Choose an example",
            list(examples_dict.keys()),
            key="example_selector"
        )

        if st.button(
            "Load Example",
            key="load_example_button"
        ):
            st.session_state.visualization_smiles = examples_dict[example]
            st.rerun()

    mol = get_molecule(st.session_state.visualization_smiles)

    if mol is None:

        st.error(
            "❌ Invalid SMILES. Please enter a valid SMILES string."
        )

    else:

        st.success("✅ Valid molecular structure detected.")

        col1, col2 = st.columns(2)

        # ----------------------------------------------------
        # 2D STRUCTURE
        # ----------------------------------------------------

        with col1:

            st.subheader("2D Molecular Structure")

            image = Draw.MolToImage(
                mol,
                size=(500, 400)
            )

            st.image(image)

        # ----------------------------------------------------
        # MOLECULAR INFORMATION
        # ----------------------------------------------------

        with col2:

            st.subheader("Molecular Information")

            properties = calculate_properties(mol)

            st.write(
                "**Canonical SMILES:**",
                Chem.MolToSmiles(mol)
            )

            st.write(
                "**Molecular Formula:**",
                properties["Molecular Formula"]
            )

            st.write(
                "**Number of atoms:**",
                mol.GetNumAtoms()
            )

            st.write(
                "**Number of bonds:**",
                mol.GetNumBonds()
            )

        st.divider()

        # ----------------------------------------------------
        # 3D VISUALIZATION
        # ----------------------------------------------------

        st.subheader("🌐 3D Molecular Visualization")

        if st.button(
            "Generate 3D Structure",
            key="generate_3d_button"
        ):

            with st.spinner("Generating 3D molecular structure..."):

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
                            height=520,
                            scrolling=False
                        )

                    else:

                        st.warning(
                            "Unable to generate a 3D structure "
                            "for this molecule."
                        )

                except Exception as e:

                    st.error(
                        f"Error generating 3D structure: {e}"
                    )


# ============================================================
# PAGE 4: MOLECULAR DESCRIPTORS
# ============================================================

elif page == "📊 Molecular Descriptors":

    st.title("📊 Molecular Descriptor Calculator")

    st.markdown("""
Calculate important physicochemical descriptors from a
molecular SMILES representation.
""")

    smiles = st.text_input(
        "Enter Molecular SMILES",
        value="CCO",
        key="descriptor_smiles"
    )

    mol = get_molecule(smiles)

    if mol is None:

        st.error("❌ Invalid SMILES.")

    else:

        properties = calculate_properties(mol)

        # ----------------------------------------------------
        # METRICS
        # ----------------------------------------------------

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

        st.divider()

        # ----------------------------------------------------
        # TABLE
        # ----------------------------------------------------

        descriptor_data = pd.DataFrame({
            "Descriptor": [
                "Molecular Weight",
                "LogP",
                "TPSA",
                "Hydrogen Bond Donors",
                "Hydrogen Bond Acceptors",
                "Rotatable Bonds",
                "Ring Count",
                "Molecular Formula"
            ],
            "Value": [
                properties["Molecular Weight"],
                properties["LogP"],
                properties["TPSA"],
                properties["HBD"],
                properties["HBA"],
                properties["Rotatable Bonds"],
                properties["Ring Count"],
                properties["Molecular Formula"]
            ]
        })

        st.subheader("Complete Results")

        st.dataframe(
            descriptor_data,
            use_container_width=True,
            hide_index=True
        )

        st.divider()

        # ----------------------------------------------------
        # INTERPRETATION
        # ----------------------------------------------------

        st.subheader("🔍 Interpretation")

        if properties["LogP"] > 3:

            st.warning(
                "Higher LogP suggests greater lipophilicity."
            )

        else:

            st.info(
                "The molecule has low to moderate lipophilicity."
            )

        if properties["TPSA"] > 90:

            st.info(
                "Higher TPSA indicates a relatively polar molecule."
            )

        if properties["HBD"] > 0:

            st.success(
                f"The molecule contains "
                f"{properties['HBD']} hydrogen bond donor(s)."
            )

        if properties["HBA"] > 0:

            st.success(
                f"The molecule contains "
                f"{properties['HBA']} hydrogen bond acceptor(s)."
            )


# ============================================================
# PAGE 5: STRUCTURE-PROPERTY ANALYSIS
# ============================================================

elif page == "📈 Structure–Property Analysis":

    st.title("📈 Structure–Property Relationship Analysis")

    st.markdown("""
Compare multiple molecules to investigate how differences in
chemical structure influence molecular properties.
""")

    default_molecules = """Ethanol,CCO
Benzene,c1ccccc1
Phenol,Oc1ccccc1
Acetic Acid,CC(=O)O
Caffeine,CN1C=NC2=C1C(=O)N(C(=O)N2C)C"""

    molecule_input = st.text_area(
        "Enter molecules as: Name,SMILES",
        value=default_molecules,
        height=200,
        key="molecule_analysis_input"
    )

    rows = []
    invalid_molecules = []

    for line in molecule_input.strip().splitlines():

        line = line.strip()

        if not line:
            continue

        if "," not in line:

            invalid_molecules.append(line)
            continue

        name, smiles = line.split(",", 1)

        mol = get_molecule(smiles)

        if mol is not None:

            properties = calculate_properties(mol)

            rows.append({
                "Molecule": name.strip(),
                "SMILES": Chem.MolToSmiles(mol),
                "MW": properties["Molecular Weight"],
                "LogP": properties["LogP"],
                "TPSA": properties["TPSA"],
                "HBD": properties["HBD"],
                "HBA": properties["HBA"],
                "Rotatable Bonds": properties["Rotatable Bonds"]
            })

        else:

            invalid_molecules.append(
                f"{name.strip()} ({smiles.strip()})"
            )

    if invalid_molecules:

        st.warning(
            "Some entries could not be processed:\n\n"
            + "\n".join(
                [f"- {item}" for item in invalid_molecules]
            )
        )

    if rows:

        df = pd.DataFrame(rows)

        st.subheader("📊 Comparison Table")

        st.dataframe(
            df,
            use_container_width=True,
            hide_index=True
        )

        st.divider()

        st.subheader("📈 Property Comparison Graph")

        properties_list = [
            "MW",
            "LogP",
            "TPSA",
            "HBD",
            "HBA",
            "Rotatable Bonds"
        ]

        col1, col2 = st.columns(2)

        with col1:

            x_axis = st.selectbox(
                "Select X-axis",
                properties_list,
                index=0,
                key="analysis_x_axis"
            )

        with col2:

            y_axis = st.selectbox(
                "Select Y-axis",
                properties_list,
                index=1,
                key="analysis_y_axis"
            )

        fig, ax = plt.subplots(figsize=(8, 5))

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
                ),
                xytext=(5, 5),
                textcoords="offset points"
            )

        ax.set_xlabel(x_axis)
        ax.set_ylabel(y_axis)
        ax.set_title(f"{x_axis} vs {y_axis}")

        ax.grid(True, alpha=0.3)

        st.pyplot(fig)

        plt.close(fig)

        st.divider()

        st.subheader("🧠 Student Observation")

        st.markdown("""
Based on your analysis, answer:

1. Which molecule has the highest molecular weight?

2. Which molecule has the highest LogP?

3. Which molecule has the highest TPSA?

4. How do functional groups influence HBD and HBA?

5. What relationship do you observe between molecular
structure and molecular properties?
""")

    else:

        st.error(
            "No valid molecules were found. "
            "Please enter molecules in the format: Name,SMILES"
        )


# ============================================================
# PAGE 6: ASSESSMENT
# ============================================================

elif page == "📝 Assessment":

    st.title("📝 Virtual Lab Assessment")

    st.markdown("""
Test your understanding of the Cheminformatics experiment.
""")

    st.divider()

    questions = [

        {
            "question": "What does SMILES represent?",

            "options": [
                "A molecular text representation",
                "A molecular weight",
                "A spectroscopy technique",
                "A laboratory instrument"
            ],

            "answer": "A molecular text representation"
        },

        {
            "question":
                "Which descriptor is associated with lipophilicity?",

            "options": [
                "TPSA",
                "LogP",
                "HBD",
                "HBA"
            ],

            "answer": "LogP"
        },

        {
            "question": "TPSA is mainly related to:",

            "options": [
                "Molecular polarity",
                "Atomic number",
                "Colour",
                "Temperature"
            ],

            "answer": "Molecular polarity"
        },

        {
            "question": "HBD stands for:",

            "options": [
                "Hydrogen Bond Donor",
                "High Bond Density",
                "Hydrogen Binary Data",
                "Heavy Bond Donor"
            ],

            "answer": "Hydrogen Bond Donor"
        },

        {
            "question":
                "Which Python toolkit is used in this virtual lab?",

            "options": [
                "RDKit",
                "MS Word",
                "PowerPoint",
                "Photoshop"
            ],

            "answer": "RDKit"
        }
    ]

    answers = []

    for i, item in enumerate(questions):

        answer = st.radio(
            f"{i + 1}. {item['question']}",
            item["options"],
            key=f"assessment_question_{i}"
        )

        answers.append(answer)

    st.divider()

    if st.button(
        "Submit Assessment",
        key="submit_assessment_button"
    ):

        score = 0

        for i, answer in enumerate(answers):

            if answer == questions[i]["answer"]:

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

        if percentage == 100:

            st.balloons()

            st.success(
                "🎉 Excellent! You successfully completed "
                "the virtual laboratory."
            )

        elif percentage >= 60:

            st.info(
                "👍 Good work! Review the experiment to improve "
                "your understanding."
            )

        else:

            st.warning(
                "Please review the theory and repeat the experiment."
            )
