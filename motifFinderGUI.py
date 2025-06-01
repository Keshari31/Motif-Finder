# Importing the required modules from the tkinter library
import tkinter as tk
from tkinter import filedialog, messagebox
from functools import partial

# Define a class for the Motif Finder GUI
class MotifFinderGUI :
    def __init__(self, root):
        # Initialize the GUI with the given root window
        self.root = root
        # Set the title of the root window
        self.root.title("Motif Finder")

        # Initialize the GUI components
        self.initialize_gui()

    # Method to initialize the GUI components
    def initialize_gui(self):
        # Frame 1: Menu Frame
        menu_frame = tk.Frame(self.root)
        menu_frame.grid(row=0, column=0, sticky="nsew", padx=0, pady=5)

        tk.Label(menu_frame, text="Select an option:").pack(pady=5)

        # tk.Button(menu_frame, text="Sequence Element", command=self.seq_element).pack(pady=5)
        # tk.Button(menu_frame, text="Motif Pattern", command=self.motif_pattern).pack(pady=5)
        # tk.Button(menu_frame, text="Predict Promoter", command=self.predict_promoter).pack(pady=5)
        tk.Button(menu_frame, text="Sequence Type", command=self.sequence_type).pack(pady=5)

        # Frame 2: Information Frame
        info_frame = tk.Frame(self.root, bg="lightgray")
        info_frame.grid(row=0, column=1, sticky="nsew", padx=10, pady=10)

        tk.Label(info_frame, text="Information Section", font=("Helvetica", 14)).pack(pady=10)

        # Frame 3: Output Frame
        output_frame = tk.Frame(self.root, bg="lightblue")
        output_frame.grid(row=2, column=1, sticky="nsew", padx=10, pady=10)

        tk.Label(output_frame, text="Output Section", font=("Helvetica", 14)).pack(pady=10)

        # Configure grid weights to make the frames expandable
        self.root.grid_rowconfigure(0, weight=1)
        self.root.grid_columnconfigure(0, weight=1)
        self.root.grid_columnconfigure(1, weight=1)
        self.root.grid_columnconfigure(2, weight=1)

    # Method to handle the 'Sequence Type' menu option
    def sequence_type(self):
        self.sequence_type = ""
        # Check the sequence against DNA, RNA, and protein alphabets
        if all(char.upper() in self.DNA_ALPHABET for char in self.sequence):
            self.sequence_type = 'DNA'
        elif all(char.upper() in self.RNA_ALPHABET for char in self.sequence):
            self.sequence_type = 'RNA'
        elif all(char.upper() in self.PROTEIN_ALPHABET for char in self.sequence):
            self.sequence_type = 'Amino acid'
        else:
            self.sequence_type = 'unknown'
        return self.sequence_type
        pass

# Initialize Tkinter and create the main window
root = tk.Tk()
root.geometry("800x500")
root.title("MOTIF FINDER")
# Create an instance of the MotifFinderGUI class with the root window
app = MotifFinderGUI(root)
# Start the Tkinter event loop to display the GUI
root.mainloop()





