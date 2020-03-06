import uproot

print(uproot.open("nano.root")["Events"].keys())
