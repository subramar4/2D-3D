import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import Crippen
from rdkit.Chem import Lipinski
from rdkit.Chem import rdMolDescriptors

import py3Dmol
import streamlit.components.v1 as components


# ============================================================
# PAGE CONFIGURATION
# ============================================================

st.set_page_config(
    page_title="Cheminformatics Virtual Laboratory",
    page_icon="🧪",
    layout="wide"
)


# ============================================================
# SIDEBAR NAVIGATION
# ============================================================

st.sidebar.title("🧪 Virtual Lab")

st.sidebar.markdown("---")

page = st.sidebar.radio(
    "Navigation",
    [
        "🏠 Home",
        "📖 Theory",
        "🧬 SMILES & Molecular Visualization",
        "📊 Molecular Descriptors",
        "📈 Structure–Property Analysis",
        "📝 Assessment"
    ]
)

st.sidebar.markdown("---")

st.sidebar.info("""
### Learning Objectives

• SMILES representation

• 2D and 3D visualization

• Molecular descriptors

• Structure–property relationships
""")


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

    This experiment introduces students to the fundamental concepts
    of **Cheminformatics** using computational tools.

    Students will learn how to represent chemical structures using
    **SMILES notation**, visualize molecules in **2D and 3D**, calculate
    important molecular descriptors, and investigate
    **structure–property relationships**.
    """)

    st.markdown("## 🎯 Learning Objectives")

    st.success("""
    ### 1. Molecular Representation

    Convert molecular representations using **SMILES notation**.
    """)

    st.success("""
    ### 2. Molecular Visualization

    Visualize chemical structures in **2D and 3D**.
    """)

    st.success("""
    ### 3. Molecular Descriptor Calculation

    Calculate important molecular descriptors:

    - Molecular Weight
    - LogP
    - TPSA
    - Hydrogen Bond Donors (HBD)
    - Hydrogen Bond Acceptors (HBA)
    """)

    st.success("""
    ### 4. Structure–Property Relationships

    Understand how molecular structure influences physical and
    chemical properties using computational methods.
    """)

    st.markdown("## 🔬 Virtual Laboratory Workflow")

    col1, col2, col3, col4 = st.columns(4)

    with col1:
        st.info("""
        ### Step 1

        📖 Learn Theory
        """)

    with col2:
        st.info("""
        ### Step 2

        🧬 Enter SMILES
        """)

    with col3:
        st.info("""
        ### Step 3

        📊 Calculate Properties
        """)

    with col4:
        st.info("""
        ### Step 4

        📈 Analyze Results
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

    st.markdown("## 🧬 Molecular Representation")

    st.markdown("""
    Molecules can be represented in several ways.

    ### Molecular Formula

    Shows the number of atoms.

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

    st.markdown("## 📊 Molecular Descriptors")

    st.markdown("""
    Molecular descriptors are numerical values that describe
    important molecular properties.

    | Descriptor | Meaning |
    |---|---|
    | Molecular Weight | Molecular mass |
    | LogP | Lipophilicity |
    | TPSA | Molecular polarity |
    | HBD | Hydrogen bond donors |
    | HBA | Hydrogen bond acceptors |
    """)

    st.markdown("## 🔬 Structure–Property Relationship")

    st.markdown("""
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

    col1, col2 = st.columns([2, 1])

    with col1:

        smiles = st.text_input(
            "Enter SMILES",
            value="CCO"
        )

    with col2:

        example = st.selectbox(
            "Choose an example",
            [
                "Ethanol",
                "Benzene",
                "Acetic Acid",
                "Caffeine"
            ]
        )

    examples_dict = {

        "Ethanol": "CCO",

        "Benzene": "c1ccccc1",

        "Acetic Acid": "CC(=O)O",

        "Caffeine":
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    }

    if st.button("Load Example"):

        smiles = examples_dict[example]

        st.info(f"Selected SMILES: {smiles}")

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:

        st.error("❌ Invalid SMILES. Please enter a valid SMILES string.")

    else:

        st.success("✅ Valid molecular structure detected.")

        col1, col2 = st.columns(2)

        # -------------------------
        # 2D STRUCTURE
        # -------------------------

        with col1:

            st.subheader("2D Molecular Structure")

            image = Draw.MolToImage(
                mol,
                size=(500, 400)
            )

            st.image(image)

        # -------------------------
        # INFORMATION
        # -------------------------

        with col2:

            st.subheader("Molecular Information")

            canonical_smiles = Chem.MolToSmiles(mol)

            st.write(
                "**Canonical SMILES:**",
                canonical_smiles
            )

            st.write(
                "**Number of atoms:**",
                mol.GetNumAtoms()
            )

            st.write(
                "**Number of bonds:**",
                mol.GetNumBonds()
            )

            st.write(
                "**Molecular formula:**",
                rdMolDescriptors.CalcMolFormula(mol)
            )

        st.divider()

        # -------------------------
        # 3D VISUALIZATION
        # -------------------------

        st.subheader("🌐 3D Molecular Visualization")

        if st.button("Generate 3D Structure"):

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

                view.setStyle(
                    {"stick": {}}
                )

                view.setBackgroundColor("white")

                view.zoomTo()

                components.html(
                    view._make_html(),
                    height=520,
                    scrolling=False
                )

            else:

                st.warning(
                    "Unable to generate a 3D structure."
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
        value="CCO"
    )

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:

        st.error("❌ Invalid SMILES.")

    else:

        # -------------------------
        # CALCULATIONS
        # -------------------------

        molecular_weight = round(
            Descriptors.MolWt(mol),
            2
        )

        logp = round(
            Crippen.MolLogP(mol),
            2
        )

        tpsa = round(
            rdMolDescriptors.CalcTPSA(mol),
            2
        )

        hbd = Lipinski.NumHDonors(mol)

        hba = Lipinski.NumHAcceptors(mol)

        # -------------------------
        # METRICS
        # -------------------------

        col1, col2, col3, col4, col5 = st.columns(5)

        col1.metric(
            "Molecular Weight",
            molecular_weight
        )

        col2.metric(
            "LogP",
            logp
        )

        col3.metric(
            "TPSA",
            tpsa
        )

        col4.metric(
            "HBD",
            hbd
        )

        col5.metric(
            "HBA",
            hba
        )

        st.divider()

        # -------------------------
        # TABLE
        # -------------------------

        descriptor_data = pd.DataFrame({

            "Descriptor": [
                "Molecular Weight",
                "LogP",
                "TPSA",
                "Hydrogen Bond Donors",
                "Hydrogen Bond Acceptors"
            ],

            "Value": [
                molecular_weight,
                logp,
                tpsa,
                hbd,
                hba
            ]
        })

        st.subheader("Complete Results")

        st.dataframe(
            descriptor_data,
            use_container_width=True,
            hide_index=True
        )

        # -------------------------
        # INTERPRETATION
        # -------------------------

        st.subheader("🔍 Interpretation")

        if logp > 3:

            st.warning(
                "Higher LogP suggests greater lipophilicity."
            )

        else:

            st.info(
                "The molecule has low to moderate lipophilicity."
            )

        if tpsa > 90:

            st.info(
                "Higher TPSA indicates a relatively polar molecule."
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
        height=200
    )

    rows = []

    for line in molecule_input.strip().splitlines():

        if "," in line:

            name, smiles = line.split(",", 1)

            mol = Chem.MolFromSmiles(
                smiles.strip()
            )

            if mol is not None:

                rows.append({

                    "Molecule": name.strip(),

                    "MW":
                    round(
                        Descriptors.MolWt(mol),
                        2
                    ),

                    "LogP":
                    round(
                        Crippen.MolLogP(mol),
                        2
                    ),

                    "TPSA":
                    round(
                        rdMolDescriptors.CalcTPSA(mol),
                        2
                    ),

                    "HBD":
                    Lipinski.NumHDonors(mol),

                    "HBA":
                    Lipinski.NumHAcceptors(mol)
                })

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

        properties = [
            "MW",
            "LogP",
            "TPSA",
            "HBD",
            "HBA"
        ]

        col1, col2 = st.columns(2)

        with col1:

            x_axis = st.selectbox(
                "Select X-axis",
                properties,
                index=0
            )

        with col2:

            y_axis = st.selectbox(
                "Select Y-axis",
                properties,
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

        st.pyplot(fig)

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


# ============================================================
# PAGE 6: ASSESSMENT
# ============================================================

elif page == "📝 Assessment":

    st.title("📝 Virtual Lab Assessment")

    st.markdown("""
    Test your understanding of the Cheminformatics experiment.
    """)

    questions = [

        {
            "question":
            "What does SMILES represent?",

            "options": [
                "A molecular text representation",
                "A molecular weight",
                "A spectroscopy technique",
                "A laboratory instrument"
            ],

            "answer":
            "A molecular text representation"
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

            "answer":
            "LogP"
        },

        {
            "question":
            "TPSA is mainly related to:",

            "options": [
                "Molecular polarity",
                "Atomic number",
                "Colour",
                "Temperature"
            ],

            "answer":
            "Molecular polarity"
        },

        {
            "question":
            "HBD stands for:",

            "options": [
                "Hydrogen Bond Donor",
                "High Bond Density",
                "Hydrogen Binary Data",
                "Heavy Bond Donor"
            ],

            "answer":
            "Hydrogen Bond Donor"
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

            "answer":
            "RDKit"
        }
    ]

    answers = []

    for i, item in enumerate(questions):

        answer = st.radio(

            f"{i+1}. {item['question']}",

            item["options"],

            key=f"question_{i}"
        )

        answers.append(answer)

    if st.button("Submit Assessment"):

        score = 0

        for i, answer in enumerate(answers):

            if answer == questions[i]["answer"]:

                score += 1

        percentage = (
            score / len(questions)
        ) * 100

        st.divider()

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
                "🎉 Excellent! You successfully completed the virtual laboratory."
            )

        elif percentage >= 60:

            st.info(
                "👍 Good work! Review the experiment to improve your understanding."
            )

        else:

            st.warning(
                "Please review the theory and repeat the experiment."
            )
```
