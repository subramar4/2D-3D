from reportlab.lib.pagesizes import A4
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY
from reportlab.lib.colors import HexColor
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak
from reportlab.lib.units import inch

path = "/mnt/data/Cheminformatics_Virtual_Laboratory_Final_Report.pdf"

styles = getSampleStyleSheet()
title = ParagraphStyle("title", parent=styles["Title"], fontSize=24, leading=30, alignment=TA_CENTER, textColor=HexColor("#17365D"))
sub = ParagraphStyle("sub", parent=styles["Normal"], fontSize=13, leading=19, alignment=TA_CENTER)
h1 = ParagraphStyle("h1", parent=styles["Heading1"], fontSize=16, leading=20, textColor=HexColor("#17365D"), spaceBefore=12, spaceAfter=8)
h2 = ParagraphStyle("h2", parent=styles["Heading2"], fontSize=12.5, leading=16, textColor=HexColor("#2F5597"), spaceBefore=8, spaceAfter=5)
body = ParagraphStyle("body", parent=styles["BodyText"], fontSize=10.5, leading=15, alignment=TA_JUSTIFY, spaceAfter=7)
bullet = ParagraphStyle("bullet", parent=body, leftIndent=18, firstLineIndent=-10, spaceAfter=4)

doc = SimpleDocTemplate(path, pagesize=A4, leftMargin=50, rightMargin=50, topMargin=55, bottomMargin=50)
story = []

def P(text, style=body):
    story.append(Paragraph(text, style))

def bullets(items):
    for x in items:
        P("• " + x, bullet)

def table(data, widths, header=True):
    t = Table(data, colWidths=widths, repeatRows=1 if header else 0)
    cmds = [
        ("GRID",(0,0),(-1,-1),0.5,colors.grey),
        ("VALIGN",(0,0),(-1,-1),"TOP"),
        ("PADDING",(0,0),(-1,-1),7),
    ]
    if header:
        cmds += [("BACKGROUND",(0,0),(-1,0),HexColor("#17365D")),
                 ("TEXTCOLOR",(0,0),(-1,0),colors.white)]
    t.setStyle(TableStyle(cmds))
    story.append(t)

# Cover
story += [Spacer(1, 1.1*inch), Paragraph("FINAL PROJECT REPORT", sub), Spacer(1, .2*inch),
          Paragraph("Cheminformatics Virtual Laboratory", title),
          Spacer(1, .2*inch),
          Paragraph("A Streamlit-Based Interactive Platform for Molecular Representation, Visualization, Descriptor Calculation and Structure–Property Analysis", sub),
          Spacer(1, .65*inch)]
table([
    ["Project Type","Interactive Virtual Laboratory"],
    ["Platform","Streamlit Community Cloud"],
    ["Programming Language","Python"],
    ["Cheminformatics Toolkit","RDKit"],
    ["Visualization","RDKit 2D and py3Dmol 3D"],
], [2.1*inch, 3.6*inch], header=False)
story += [Spacer(1,.8*inch), Paragraph("Final implementation and deployment report", sub), PageBreak()]

P("1. Abstract", h1)
P("The Cheminformatics Virtual Laboratory was developed as an interactive web-based learning environment for introducing fundamental concepts of computational chemistry and cheminformatics. The application enables users to study theoretical concepts, enter molecular structures using SMILES notation, visualize molecules in two and three dimensions, calculate important molecular descriptors, compare structure–property relationships, and complete an assessment. The complete application was implemented in a single Python file (app.py) using Streamlit for the user interface and RDKit for chemical structure processing. Pandas, Matplotlib and py3Dmol support data handling, graphical analysis and three-dimensional visualization. The application was deployed on Streamlit Community Cloud after resolving Python and Linux dependency issues.")

P("2. Objectives", h1)
bullets([
    "Develop an interactive virtual laboratory for teaching basic cheminformatics concepts.",
    "Introduce molecular representation using SMILES notation.",
    "Generate and display two-dimensional and three-dimensional molecular structures.",
    "Calculate important physicochemical molecular descriptors.",
    "Investigate structure–property relationships using molecular data.",
    "Assess student understanding through an integrated quiz module.",
    "Deploy the completed application as a web-accessible Streamlit application."
])

P("3. Software and Technologies Used", h1)
table([
    ["Technology","Purpose"],
    ["Python","Core programming language."],
    ["Streamlit","Interactive web interface and application deployment."],
    ["RDKit","SMILES processing, structure generation and descriptor calculation."],
    ["Pandas","Tabular organization and presentation of molecular data."],
    ["Matplotlib","Structure–property comparison graphs."],
    ["py3Dmol","Interactive three-dimensional molecular visualization."],
    ["Streamlit Community Cloud","Cloud deployment platform."]
], [1.65*inch,4.05*inch])

P("4. System Design and Architecture", h1)
P("The project follows a single-file Streamlit architecture. All major modules are contained in app.py, while sidebar navigation controls which component is displayed. Only one st.set_page_config() call is used, preventing page-configuration conflicts in a single-file application.")
table([
    ["USER"],
    ["STREAMLIT WEB INTERFACE"],
    ["SIDEBAR NAVIGATION"],
    ["Home | Theory | Visualization | Descriptors | Analysis | Assessment"],
    ["RDKit | Pandas | Matplotlib | py3Dmol"],
    ["Interactive Learning Results"]
], [5.7*inch])
story.append(PageBreak())

P("5. Functional Modules", h1)
modules = [
("5.1 Home Module","Introduces the virtual laboratory, learning objectives, workflow and example molecules."),
("5.2 Theory Module","Provides introductory material on cheminformatics, molecular representation, SMILES notation, molecular descriptors and structure–property relationships."),
("5.3 SMILES and Molecular Visualization","Allows users to enter or select a SMILES string, validate the molecular structure and generate a two-dimensional molecular image."),
("5.4 Three-Dimensional Visualization","Uses RDKit to prepare a molecular geometry and py3Dmol to display the molecule interactively in three dimensions."),
("5.5 Molecular Descriptor Calculator","Calculates Molecular Weight, LogP, TPSA, HBD, HBA, Rotatable Bonds, Ring Count and Molecular Formula."),
("5.6 Structure–Property Analysis","Compares multiple molecules in a table and generates selectable scatter plots for molecular property analysis."),
("5.7 Assessment Module","Provides multiple-choice questions, calculates the score and displays feedback to the learner.")
]
for a,b in modules:
    P(a,h2); P(b)

P("6. Implementation Details", h1)
P("The application uses a sidebar radio navigation system to select each module. RDKit converts valid SMILES strings into molecular objects. The program checks molecular validity before visualization and property calculations. Pandas organizes calculated results, while Matplotlib provides graphical comparison of selected molecular properties.")

P("6.1 Molecular Descriptors", h2)
table([
    ["Descriptor","Interpretation"],
    ["Molecular Weight","Total molecular mass."],
    ["LogP","Descriptor associated with molecular lipophilicity."],
    ["TPSA","Topological polar surface area associated with polarity."],
    ["HBD","Number of hydrogen bond donor sites."],
    ["HBA","Number of hydrogen bond acceptor sites."],
    ["Rotatable Bonds","Indicator of molecular flexibility."],
    ["Ring Count","Number of rings identified in the molecular structure."]
], [1.7*inch,4.0*inch])

P("7. Deployment and Dependency Configuration", h1)
P("The completed application was deployed using Streamlit Community Cloud. During deployment, dependency-related errors were identified through application logs and corrected by configuring the required Python packages and Linux system libraries.")

P("7.1 requirements.txt", h2)
table([["Package"],["streamlit"],["pandas"],["matplotlib"],["py3Dmol"],["rdkit"]],[3.0*inch])

P("7.2 packages.txt", h2)
P("The RDKit drawing module required Linux libraries for graphical rendering. The missing libXrender.so.1 error was addressed through packages.txt.")
table([["System Package"],["libxrender1"],["libxext6"],["libsm6"]],[3.0*inch])

story.append(PageBreak())
P("8. Problems Encountered and Solutions", h1)
table([
    ["Problem","Cause","Solution"],
    ["Matplotlib module error","Required dependency unavailable.","Added Matplotlib to requirements.txt."],
    ["RDKit import error","Required supporting libraries unavailable.","Configured the required dependencies."],
    ["libXrender.so.1 missing","Linux X rendering library was absent.","Added libxrender1, libxext6 and libsm6 through packages.txt."],
    ["Single-file configuration conflict","Multiple page configurations can conflict.","Used one st.set_page_config() call and sidebar routing."]
], [1.55*inch,2.0*inch,2.15*inch])

P("9. Final Outcome", h1)
P("The Cheminformatics Virtual Laboratory was successfully developed as a functional interactive web application. The final system integrates theoretical learning, molecular representation, two-dimensional visualization, three-dimensional visualization, molecular descriptor calculation, structure–property comparison and student assessment within one Streamlit application.")
bullets([
    "Single-file app.py implementation for the complete application.",
    "SMILES input and molecular validation.",
    "RDKit-based two-dimensional molecular visualization.",
    "Three-dimensional visualization using RDKit and py3Dmol.",
    "Automatic calculation of multiple physicochemical descriptors.",
    "Graphical comparison of molecular properties.",
    "Integrated assessment with automatic scoring.",
    "Successful deployment configuration after resolving dependency issues."
])

P("10. Educational Significance", h1)
P("The virtual laboratory provides a practical environment for learning computational chemistry concepts. Students can directly observe how a textual chemical representation is converted into a molecular structure and how structural differences influence calculated molecular properties. Integrating theory, visualization, calculation, analysis and assessment within one platform supports interactive and self-directed learning.")

P("11. Future Scope", h1)
bullets([
    "Add molecular file upload support for MOL and SDF formats.",
    "Expand the molecular descriptor library.",
    "Provide downloadable calculation results and laboratory reports.",
    "Include additional structure–property prediction models.",
    "Expand the assessment question bank.",
    "Add student progress tracking.",
    "Introduce additional molecular visualization and analysis tools."
])

P("12. Conclusion", h1)
P("The Cheminformatics Virtual Laboratory demonstrates the successful integration of Python-based scientific computing and web application development for chemical education. Streamlit and RDKit were combined to create an interactive learning environment for molecular representation, visualization, descriptor calculation and structure–property analysis. The deployment issues encountered during development were systematically resolved through dependency configuration, resulting in a functional application suitable for educational use.")

def footer(canvas, doc):
    canvas.saveState()
    canvas.setFont("Helvetica", 8.5)
    canvas.setFillColor(colors.grey)
    canvas.drawCentredString(A4[0]/2, 28, f"Cheminformatics Virtual Laboratory – Final Project Report | Page {doc.page}")
    canvas.restoreState()

doc.build(story, onFirstPage=footer, onLaterPages=footer)
print(path)
