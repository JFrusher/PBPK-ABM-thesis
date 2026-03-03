import xml.etree.ElementTree as ET

tree = ET.parse('PhysiCell_settings.xml')
root = tree.getroot()

cds = root.findall('cell_definitions//cell_definition//custom_data')
# print("cds=",cds)
for cd_elm in cds:
    print(f"Parent Tag: {cd_elm.tag}, Attributes: {cd_elm.attrib}")
    # Iterate through the children of the current 'parent' element
    num_vars = 0
    for child in cd_elm:
        print(f"  Child Tag: {child.tag}, Text: {child.text}")
        num_vars += 1
    print("--- num_vars=",num_vars)
