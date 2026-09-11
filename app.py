import streamlit as st

# -------------------------------------------------
# PAGE CONFIGURATION
# -------------------------------------------------
st.set_page_config(
    page_title="Cheminformatics Virtual Lab",
    page_icon="🧪",
    layout="wide"
)

# -------------------------------------------------
# HEADER
# -------------------------------------------------
st.title("🧪 Cheminformatics Virtual Laboratory")

st.subheader(
    "Computational Representation, Visualization and Analysis "
    "of Chemical Structures"
)

st.divider()

# -------------------------------------------------
# INTRODUCTION
# -------------------------------------------------
st.markdown("""
## Welcome to the Virtual Laboratory

This virtual laboratory introduces students to the fundamental concepts
of **Cheminformatics** using computational tools.

Students can represent chemical structures using **SMILES notation**,
visualize molecules in **2D and 3D**, calculate important molecular
descriptors, and investigate **structure–property relationships**.
""")

# -------------------------------------------------
# LEARNING OBJECTIVES
# -------------------------------------------------
st.markdown("## 🎯 Learning Objectives")

st.markdown("""
By completing this virtual experiment, students will be able to:

✅ Convert molecular representations using **SMILES notation**.

✅ Visualize chemical structures in **2D and 3D**.

✅ Calculate important molecular descriptors:

- Molecular Weight
- LogP
- TPSA
- Hydrogen Bond Donors (HBD)
- Hydrogen Bond Acceptors (HBA)

✅ Understand **structure–property relationships** using
computational methods.
""")

# -------------------------------------------------
# HOW TO USE
# -------------------------------------------------
st.markdown("## 🧭 Virtual Lab Workflow")

col1, col2, col3, col4 = st.columns(4)

with col1:
    st.info("""
    ### 1️⃣ Theory

    Learn the basics of Cheminformatics and molecular representations.
    """)

with col2:
    st.success("""
    ### 2️⃣ Visualization

    Convert SMILES into 2D and 3D molecular structures.
    """)

with col3:
    st.warning("""
    ### 3️⃣ Analysis

    Calculate molecular descriptors and properties.
    """)

with col4:
    st.error("""
    ### 4️⃣ Assessment

    Test your understanding using the quiz.
    """)

# -------------------------------------------------
# EXAMPLE MOLECULES
# -------------------------------------------------
st.markdown("## 🧬 Example SMILES")

examples = {
    "Ethanol": "CCO",
    "Benzene": "c1ccccc1",
    "Acetic Acid": "CC(=O)O",
    "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
}

for name, smiles in examples.items():
    st.code(f"{name}: {smiles}")

st.divider()

st.success(
    "Use the navigation menu on the left to begin the Cheminformatics Virtual Laboratory."
)
import streamlit as st

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

import py3Dmol
import streamlit.components.v1 as components


st.set_page_config(
    page_title="Molecular Visualization",
    page_icon="🧪",
    layout="wide"
)

st.title("🧪 Molecular Representation and Visualization")

st.markdown("""
Enter a **SMILES string** to generate and visualize the molecular
structure.
""")

# -------------------------------------------------
# USER INPUT
# -------------------------------------------------

smiles = st.text_input(
    "Enter SMILES",
    value="CCO"
)

# -------------------------------------------------
# CREATE MOLECULE
# -------------------------------------------------

mol = Chem.MolFromSmiles(smiles)

if mol is None:

    st.error("❌ Invalid SMILES. Please enter a valid SMILES string.")

else:

    col1, col2 = st.columns(2)

    # ---------------------------------------------
    # 2D STRUCTURE
    # ---------------------------------------------

    with col1:

        st.subheader("2D Molecular Structure")

        image = Draw.MolToImage(
            mol,
            size=(500, 400)
        )

        st.image(image)

    # ---------------------------------------------
    # MOLECULAR INFORMATION
    # ---------------------------------------------

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

        st.success(
            "SMILES successfully converted into a molecular structure."
        )

    st.divider()

    # ---------------------------------------------
    # 3D STRUCTURE
    # ---------------------------------------------

    st.subheader("🧬 3D Molecular Visualization")

    if st.button("Generate 3D Structure"):

        mol3d = Chem.AddHs(mol)

        status = AllChem.EmbedMolecule(
            mol3d,
            randomSeed=42
        )

        if status == 0:

            try:

                AllChem.MMFFOptimizeMolecule(mol3d)

            except:

                pass

            mol_block = Chem.MolToMolBlock(mol3d)

            view = py3Dmol.view(
                width=800,
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
                "Unable to generate the 3D structure for this molecule."
            )
            import streamlit as st
import pandas as pd

from rdkit import Chem
from rdkit.Chem import (
    Descriptors,
    Crippen,
    Lipinski,
    rdMolDescriptors
)


st.set_page_config(
    page_title="Molecular Descriptors",
    page_icon="📊",
    layout="wide"
)

st.title("📊 Molecular Descriptor Calculator")

st.markdown("""
Enter a molecular SMILES string to calculate important
physicochemical descriptors.
""")

# -------------------------------------------------
# INPUT
# -------------------------------------------------

smiles = st.text_input(
    "Enter SMILES",
    value="CCO"
)

mol = Chem.MolFromSmiles(smiles)

if mol is None:

    st.error("❌ Invalid SMILES.")

else:

    # ---------------------------------------------
    # CALCULATE DESCRIPTORS
    # ---------------------------------------------

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

    rotatable_bonds = Lipinski.NumRotatableBonds(mol)

    ring_count = Lipinski.RingCount(mol)

    # ---------------------------------------------
    # METRICS
    # ---------------------------------------------

    col1, col2, col3 = st.columns(3)

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

    st.divider()

    # ---------------------------------------------
    # TABLE
    # ---------------------------------------------

    data = {

        "Descriptor": [
            "Molecular Weight",
            "LogP",
            "TPSA",
            "Hydrogen Bond Donors",
            "Hydrogen Bond Acceptors",
            "Rotatable Bonds",
            "Ring Count"
        ],

        "Value": [
            molecular_weight,
            logp,
            tpsa,
            hbd,
            hba,
            rotatable_bonds,
            ring_count
        ]
    }

    df = pd.DataFrame(data)

    st.subheader("Complete Molecular Descriptor Table")

    st.dataframe(
        df,
        use_container_width=True,
        hide_index=True
    )

    # ---------------------------------------------
    # INTERPRETATION
    # ---------------------------------------------

    st.subheader("🔍 Interpretation")

    if logp > 3:

        st.warning(
            "The molecule has relatively high lipophilicity."
        )

    else:

        st.info(
            "The molecule has low to moderate lipophilicity."
        )

    if tpsa > 90:

        st.info(
            "The molecule has a relatively high polar surface area."
        )

    if hbd > 0:

        st.success(
            f"The molecule contains {hbd} hydrogen bond donor(s)."
        )

    if hba > 0:

        st.success(
            f"The molecule contains {hba} hydrogen bond acceptor(s)."
        )
        import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt

from rdkit import Chem
from rdkit.Chem import (
    Descriptors,
    Crippen,
    Lipinski,
    rdMolDescriptors
)


st.set_page_config(
    page_title="Structure Property Analysis",
    page_icon="📈",
    layout="wide"
)

st.title("📈 Structure–Property Relationship Analysis")

st.markdown("""
Compare different molecules and investigate how molecular structure
influences physicochemical properties.
""")

# -------------------------------------------------
# DEFAULT MOLECULES
# -------------------------------------------------

default_data = """Ethanol,CCO
Phenol,Oc1ccccc1
Acetic Acid,CC(=O)O
Caffeine,CN1C=NC2=C1C(=O)N(C(=O)N2C)C
Benzene,c1ccccc1
"""

# -------------------------------------------------
# USER INPUT
# -------------------------------------------------

text = st.text_area(
    "Enter molecules in the format: Name,SMILES",
    value=default_data,
    height=200
)

# -------------------------------------------------
# PROCESS MOLECULES
# -------------------------------------------------

rows = []

for line in text.strip().splitlines():

    if "," not in line:
        continue

    name, smiles = line.split(",", 1)

    mol = Chem.MolFromSmiles(
        smiles.strip()
    )

    if mol is not None:

        rows.append({

            "Molecule": name.strip(),

            "SMILES":
            Chem.MolToSmiles(mol),

            "MW":
            round(Descriptors.MolWt(mol), 2),

            "LogP":
            round(Crippen.MolLogP(mol), 2),

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


# -------------------------------------------------
# DISPLAY RESULTS
# -------------------------------------------------

if rows:

    df = pd.DataFrame(rows)

    st.subheader("📊 Molecular Property Comparison")

    st.dataframe(
        df,
        use_container_width=True,
        hide_index=True
    )

    st.divider()

    # ---------------------------------------------
    # GRAPH SELECTION
    # ---------------------------------------------

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

    # ---------------------------------------------
    # PLOT
    # ---------------------------------------------

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

    st.pyplot(fig)

    # ---------------------------------------------
    # STUDENT QUESTIONS
    # ---------------------------------------------

    st.subheader("🧠 Student Analysis")

    st.markdown("""
Answer the following questions based on your results:

1. Which molecule has the highest molecular weight?

2. Which molecule has the highest LogP?

3. Which molecule has the highest TPSA?

4. How do functional groups influence HBD and HBA?

5. What relationship can you observe between molecular structure
and molecular properties?
""")

else:

    st.warning(
        "Please enter at least one valid molecule."
    )
    import streamlit as st


st.set_page_config(
    page_title="Assessment",
    page_icon="📝",
    layout="wide"
)

st.title("📝 Cheminformatics Virtual Lab Assessment")

st.markdown("""
Answer the following questions to assess your understanding
of the experiment.
""")

# -------------------------------------------------
# QUESTIONS
# -------------------------------------------------

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
        "Which descriptor is commonly associated with lipophilicity?",

        "options": [
            "TPSA",
            "LogP",
            "HBD",
            "Molecular Formula"
        ],

        "answer":
        "LogP"
    },

    {
        "question":
        "TPSA mainly provides information about:",

        "options": [
            "Molecular polarity",
            "Atomic number",
            "Colour",
            "Melting apparatus"
        ],

        "answer":
        "Molecular polarity"
    },

    {
        "question":
        "Which toolkit is commonly used for cheminformatics in Python?",

        "options": [
            "RDKit",
            "Microsoft Word",
            "PowerPoint",
            "Excel only"
        ],

        "answer":
        "RDKit"
    },

    {
        "question":
        "HBD stands for:",

        "options": [
            "Hydrogen Bond Donor",
            "Hydrogen Bond Density",
            "High Bond Distance",
            "Hydrogen Binary Data"
        ],

        "answer":
        "Hydrogen Bond Donor"
    }
]


# -------------------------------------------------
# STUDENT ANSWERS
# -------------------------------------------------

student_answers = []

for i, item in enumerate(questions):

    answer = st.radio(

        f"{i+1}. {item['question']}",

        item["options"],

        key=f"question_{i}"
    )

    student_answers.append(answer)


# -------------------------------------------------
# SUBMIT
# -------------------------------------------------

if st.button("Submit Assessment"):

    score = 0

    for i, answer in enumerate(student_answers):

        if answer == questions[i]["answer"]:

            score += 1

    st.divider()

    st.success(
        f"Your Score: {score} / {len(questions)}"
    )

    percentage = (
        score / len(questions)
    ) * 100

    st.metric(
        "Percentage",
        f"{percentage:.1f}%"
    )

    if score == len(questions):

        st.balloons()

        st.success(
            "🎉 Excellent! You have successfully completed the Cheminformatics Virtual Laboratory."
        )

    elif percentage >= 60:

        st.info(
            "Good performance! Review the theory section to strengthen your understanding."
        )

    else:

        st.warning(
            "Please review the theory and repeat the virtual experiment."
        )