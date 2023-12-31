# Load DNA sequences
from Bio import SeqIO
from sklearn.tree import DecisionTreeClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score

reference_file = "reference.fna"
sample_file = "gene3.fna"

reference_sequence = str(SeqIO.read(reference_file, "fasta").seq)
sample_sequence = str(SeqIO.read(sample_file, "fasta").seq)

# Encode DNA bases as features
features = []

for base in reference_sequence:
    if base == 'A':
        features.append([1, 0, 0, 0])
    elif base == 'C':
        features.append([0, 1, 0, 0])
    elif base == 'G':
        features.append([0, 0, 1, 0])
    elif base == 'T':
        features.append([0, 0, 0, 1])

# Create binary labels indicating variant positions
labels = [1 if ref != sample else 0 for ref, sample in zip(reference_sequence, sample_sequence)]

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(features, labels, test_size=0.2, random_state=42)

# Train the model
model = DecisionTreeClassifier()
model.fit(X_train, y_train)

# Evaluate the model
y_pred = model.predict(X_test)
accuracy = accuracy_score(y_test, y_pred)
print(f"Model Accuracy: {accuracy}")
