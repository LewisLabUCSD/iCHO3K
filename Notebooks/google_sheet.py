
#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Module Name: google_sheet.py
Author: Pablo Di Giusto
Date Created: 2023-04-06
Last Updated: 2023-05-16
Version: 1.0.0
Description: This module facilitates the integration of the Google Sheets with python by generating creating a function to
retrieve the information from a Google Sheet as a pandas dataframe and another function to update the Google Sheet from a
pandas dataframe.

"""


import json
import pandas as pd
from google.oauth2 import service_account
from googleapiclient.discovery import build
from openpyxl.utils import get_column_letter


class GoogleSheet:

    def __init__(self, sheet_id, credentials_path):
        self.sheet_id = sheet_id
        self.credentials = service_account.Credentials.from_service_account_file(
            credentials_path,
            scopes=['https://www.googleapis.com/auth/spreadsheets']
        )
        self.sheets_api = build('sheets', 'v4', credentials=self.credentials)


    def read_google_sheet(self, sheet_name):
        '''
        Function to read Google Sheet data into a DataFrame
        '''

        # Get the sheet metadata to find the range
        sheet_metadata = self.sheets_api.spreadsheets().get(
            spreadsheetId=self.sheet_id).execute()
        sheets = sheet_metadata.get('sheets', '')
        sheet = None

        for s in sheets:
            if s['properties']['title'] == sheet_name:
                sheet = s
                break

        if sheet is None:
            raise Exception(f"Sheet '{sheet_name}' not found")

        # Calculate the row count (including empty rows)
        row_count = sheet['properties']['gridProperties']['rowCount']
        column_count = sheet['properties']['gridProperties']['columnCount']

        sheet_range = f"{sheet_name}!A1:{column_count}{row_count}"

        result = self.sheets_api.spreadsheets().values().get(
            spreadsheetId=self.sheet_id, range=sheet_range).execute()
        data = result.get('values', [])
        df = pd.DataFrame(data[1:], columns=data[0])
        return df


    def update_google_sheet(self, sheet_name, df):
        # Replace NaNs with empty strings
        df = df.fillna('')

        data = [df.columns.tolist()] + df.values.tolist()

        num_rows = df.shape[0] + 1
        num_columns = df.shape[1]
        last_column = get_column_letter(num_columns)
        sheet_range = f'{sheet_name}!A1:{last_column}{num_rows}'

        # Get the current row count of the Google Sheet
        current_row_count = self.get_row_count(sheet_name)

        # Clear the specified range
        self.sheets_api.spreadsheets().values().clear(
            spreadsheetId=self.sheet_id,
            range=sheet_range).execute()

        # Update the Google Sheet with the new data
        body = {
            'range': sheet_range,
            'values': data,
            'majorDimension': 'ROWS'
        }
        self.sheets_api.spreadsheets().values().update(
            spreadsheetId=self.sheet_id,
            range=sheet_range,
            valueInputOption='RAW',
            body=body).execute()

        # Delete extra rows from the Google Sheet if needed
        if num_rows < current_row_count:
            self.delete_extra_rows(sheet_name, num_rows, current_row_count)


    def get_row_count(self, sheet_name):
        sheet_metadata = self.sheets_api.spreadsheets().get(spreadsheetId=self.sheet_id).execute()
        sheets = sheet_metadata.get('sheets', '')
        sheet = None

        for s in sheets:
            if s['properties']['title'] == sheet_name:
                sheet = s
                break

        if sheet is None:
            raise Exception(f"Sheet '{sheet_name}' not found")

        row_count = sheet['properties']['gridProperties']['rowCount']
        return row_count


    def delete_extra_rows(self, sheet_name, start_row, end_row):
        sheet_id = None
        sheet_metadata = self.sheets_api.spreadsheets().get(spreadsheetId=self.sheet_id).execute()
        sheets = sheet_metadata.get('sheets', '')

        for s in sheets:
            if s['properties']['title'] == sheet_name:
                sheet_id = s['properties']['sheetId']
                break

        if sheet_id is None:
            raise Exception(f"Sheet '{sheet_name}' not found")

        delete_range = {
            "sheetId": sheet_id,
            "dimension": "ROWS",
            "startIndex": start_row,
            "endIndex": end_row
        }

        body = {
            "requests": [
                {
                    "deleteDimension": {
                        "range": delete_range
                    }
                }
            ]
        }

        self.sheets_api.spreadsheets().batchUpdate(
            spreadsheetId=self.sheet_id,
            body=body).execute()
