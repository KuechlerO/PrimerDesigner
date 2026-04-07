from django.shortcuts import render, redirect
from django.http import HttpResponse


def documentation(request):
    return render(request, 'primer_designer_app/documentation.html')
