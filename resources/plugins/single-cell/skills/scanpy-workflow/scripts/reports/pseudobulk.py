# JSON contract for sc_py_pseudobulk.
_sci_ctx = {}
if counts_output and pb_info and pb_info.get("design"):
    _sci_ctx[str(counts_output)] = {
        "valueSemantics": "raw_counts",
        "sampleIds": [str(entry["column"]) for entry in pb_info["design"]],
        "donorBySample": {
            str(entry["column"]): str(entry["sample"]) for entry in pb_info["design"]
        },
    }
emit({"input": tool_input, **pb_info, **({"scientificContext": _sci_ctx} if _sci_ctx else {})})