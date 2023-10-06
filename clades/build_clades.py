#!/usr/bin/env python3

import pandas as pd
from io import StringIO

# https://stackoverflow.com/questions/29481088/how-can-i-tell-if-a-string-repeats-itself-in-python
def principal_period(s):
    i = (s + s).find(s, 1, -1)
    return None if i == -1 else s[:i]


def filter_telo_repeats(l):
    r_simple_repeats = []
    # for string in the list
    for s in l:
        if isinstance(s, str):
            if principal_period(s) is None:
                if s != len(s) * s[0]:
                    r_simple_repeats.append(s)

    return r_simple_repeats


# https://stackoverflow.com/questions/1630320/what-is-the-pythonic-way-to-detect-the-last-element-in-a-for-loop
def lookahead(iterable):
    """Pass through all values from the given iterable, augmented by the
    information if there are more values to come after the current one
    (True), or if it is the last value (False).
    """
    # Get an iterator and pull the first value.
    it = iter(iterable)
    try:
        last = next(it)
    except StopIteration:
        return
    # Run the iterator to exhaustion (starting from the second value).
    for val in it:
        # Report the *previous* value (more to come).
        yield last, True
        last = val
    # Report the last value.
    yield last, False

df = pd.read_table("./curated.csv", delimiter=",")

# group by order
gb_order = df.groupby(["Order"])["Telomeric repeat"].apply(list)

# now let's print some Rust!

print(
    """
    /// All the clades for which we have data.
    pub static CLADES: &[&str] = &["""
)

for i, v in gb_order.items():
    omit_pure_repeats = filter_telo_repeats(v)
    if omit_pure_repeats:
        print('        "', i, '",', sep="")

print(
    """    ];
"""
)

print(
    """
    /// A function to get a telomeric repeat sequence
    /// given a clade name.
    pub fn return_telomere_sequence(clade: &str) -> TelomereSeq {
        let result = match clade {"""
)

for i, v in gb_order.items():
    omit_pure_repeats = filter_telo_repeats(v)
    omit_pure_repeats_list = list(set(omit_pure_repeats))

    if omit_pure_repeats_list:
        omit_pure_repeats_list_string = ""

        for telomere, has_more in lookahead(omit_pure_repeats_list):
            if has_more:
                omit_pure_repeats_list_string += '"' + telomere + '",\n\t\t\t'
            else:
                omit_pure_repeats_list_string += '"' + telomere + '"'

        print(
            f"""
        \"{i}\" => TelomereSeq {{
            clade: \"{i.capitalize()}\",
            seq: Seq(Box::new(&[{omit_pure_repeats_list_string}])),
            length: {len(omit_pure_repeats_list)},
        }},
        """
        )

print(
    """
            _ => panic!("{} is not yet accounted for in this pipeline.", clade),
        };
        result
    }"""
)
