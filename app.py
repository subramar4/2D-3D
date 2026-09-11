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
