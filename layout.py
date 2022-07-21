import plotly.graph_objs as go

layout = go.Layout(
    title="Sine and cos",
    xaxis=dict(
        title='angle',
        showgrid=True,
        zeroline=True,
        showline=True,
        showticklabels=True,
        gridwidth=1
    ),
    yaxis=dict(
        showgrid=True,
        zeroline=True,
        showline=True,
        gridcolor='#bdbdbd',
        gridwidth=2,
        zerolinecolor='#969696',
        zerolinewidth=2,
        linecolor='#636363',
        linewidth=2,
        title='VALUE',
        titlefont=dict(
            family='Arial, sans-serif',
            size=18,
            color='lightgrey'
        ),
        showticklabels=True,
        tickangle=45,
        tickfont=dict(
            family='Old Standard TT, serif',
            size=14,
            color='black'
        ),
        tickmode='linear',
        tick0=0.0,
        dtick=0.25
    )
)
