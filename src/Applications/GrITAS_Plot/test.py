#!/usr/bin/env python
import common


# Regions to focus on
regTypes = ['glo','nhe','she','tro','nam']

# Instantiate Region class members
regions=[common.Region(r) for r in regTypes]

# Bundle all regions focused on into a super set
Universe=common.Universe(regions)


# Output serialization
out=open('regions.yaml','w')
myyam=common.myYML(out)


# Form the universe aspect of yaml
Universe.serialize(myyam)
