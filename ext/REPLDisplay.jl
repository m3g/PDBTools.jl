module REPLDisplay
    using PDBTools: PDBTools, StructureView, open_browser
    import REPL

    # In the REPL, which cannot render HTML, structure views are opened in the browser.
    # Environments that render HTML (VS Code, Jupyter, Pluto) use their own displays.
    function Base.display(d::REPL.REPLDisplay, v::StructureView)
        file = open_browser(v)
        io = REPL.outstream(d.repl)
        show(io, MIME"text/plain"(), v)
        print(io, "\n    Opened in the browser: $file\n")
        return nothing
    end
end
